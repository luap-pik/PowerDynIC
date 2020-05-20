"""
assuming that you load PowerDynamics beforehand
"""

## patching PD.jl with ControlledPowerGrid

using LightGraphs: AbstractGraph

# stack dynamics with higher level controller
struct ControlledPowerGrid
    graph:: G where G <: AbstractGraph
    nodes
    lines
    controller # oop function
end

# Constructor
function ControlledPowerGrid(control, pg::PowerGrid, args...)
    ControlledPowerGrid(pg.graph, pg.nodes, pg.lines, (u, p, t) -> control(u, p, t, args...))
end

import PowerDynamics: PowerGrid

PowerGrid(pg::PowerGrid) = pg
PowerGrid(pg::ControlledPowerGrid) = PowerGrid(pg.graph, pg.nodes, pg.lines)

import PowerDynamics: State, PowerGridSolution

State(pg::ControlledPowerGrid, vec) = State(PowerGrid(pg), vec)
PowerGridSolution(vec, pg::ControlledPowerGrid) = PowerGridSolution(vec, PowerGrid(pg))

import PowerDynamics: systemsize

systemsize(cpg::ControlledPowerGrid) = systemsize(PowerGrid(cpg.graph, cpg.nodes, cpg.lines))

## define ControlledPowerGrid rhs

using OrdinaryDiffEq: ODEFunction
import PowerDynamics: rhs, get_current

function rhs(cpg::ControlledPowerGrid)
    open_loop_dyn = PowerGrid(cpg.graph, cpg.nodes, cpg.lines) |> rhs
    function cpg_rhs!(du, u, p, t)
        # Get controller state
        control_output = cpg.controller(u, p, t)
        # Calculate the derivatives, passing the controller state through
        open_loop_dyn(du, u, control_output, t)
    end
    ODEFunction{true}(cpg_rhs!, mass_matrix=open_loop_dyn.mass_matrix, syms=open_loop_dyn.syms)
end

using NetworkDynamics: network_dynamics, StaticEdgeFunction

function rhs(pg::PowerGrid)
    network_dynamics(map(construct_vertex, pg.nodes), map(construct_edge, pg.lines), pg.graph; parallel=true)
end

get_current(state::State, n) = begin
    vertices = map(construct_vertex, state.grid.nodes)
    edges = map(construct_edge, state.grid.lines)
    sef = StaticEdgeFunction(vertices,edges,state.grid.graph, parallel=true)
    (e_s, e_d) = sef(state.vec, Nothing, 0)
    total_current(e_s[n], e_d[n])
end

## operation point

import PowerDynamics: find_operationpoint
using NLsolve: nlsolve, converged
using ForwardDiff: jacobian
using LinearAlgebra: eigvals
using PowerDynamics: AbstractNode
using OrdinaryDiffEq: ODEFunction
using DifferentialEquations: DynamicSS
using NLsolve: nlsolve, converged
using LinearAlgebra: pinv

struct RootRhs_ic
    rhs
    mpm
end
function (rr::RootRhs_ic)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.)
    rr.mpm * dx .- dx
end

function RootRhs_ic(of::ODEFunction)
    mm = of.mass_matrix
    @assert mm != nothing
    mpm = pinv(mm) * mm
    RootRhs_ic(of.f, mpm)
end

"""
Makes an type specific initial guess to help the operation point search.
The voltage is of all nodes is fixed to the voltage of the first SlackAlgebraic
in the system. The other symbols are set to zero.
Inputs:
    pg: Power grid, a graph containg nodes and lines
Outputs:
    guess: Type specific initial guess
"""
function initial_guess(pg)
    if SlackAlgebraic ∉ pg.nodes .|> typeof
        @warn "There is no slack bus in the system to balance powers."
    end

    sl = findfirst(SlackAlgebraic  ∈  pg.nodes .|> typeof)
    slack = pg.nodes[sl]

    type_guesses = guess.(pg.nodes, slack.U)
    return vcat(type_guesses...)
end

function initial_guess(pg, v_vector)
    type_guesses = guess.(pg.nodes, v_vector)
    return vcat(type_guesses...)
end


"""
guess(::Type{SlackAlgebraic}) = [slack.U, 0.]         #[:u_r, :u_i]
guess(::Type{PQAlgebraic}) = [slack.U, 0.]            #[:u_r, :u_i]
guess(::Type{ThirdOrderEq}) = [slack.U, 0., 0., 0.]   #[:u_r, :u_i, :θ, :ω]
"""
function guess(n::AbstractNode, slack_voltage)
    voltage = zeros(dimension(n))
    voltage[1] = real(slack_voltage)  #[:u_r, :u_i]
    voltage[2] = imag(slack_voltage)
    return voltage
end

function guess(n::SlackAlgebraic, slack_voltage)
    @assert n.U ≈ slack_voltage
    return [real(n.U), imag(n.U)]
end

function find_steady_state(pg, op_guess=nothing)
    if op_guess isa Nothing
        op_guess = initial_guess(pg)
    end
    ode = rhs(pg)
    op_prob = ODEProblem(ode, op_guess, (0.0, 1.0))
    sol = solve(
        SteadyStateProblem(op_prob),
        DynamicSS(Rodas5(); abstol = 1e-8, reltol = 1e-6, tspan = Inf),
    )
    return State(pg, sol.u)
end

function find_operationpoint(pg::ControlledPowerGrid, ic_guess = nothing)
    if SlackAlgebraic ∉ pg.nodes .|> typeof
        @warn "There is no slack bus in the system to balance powers. Currently not making any checks concerning assumptions of whether its possible to find the fixed point"
    end
    if SwingEq ∈ pg.nodes .|> typeof
        throw(OperationPointError("found SwingEq node but these should be SwingEqLVS (just SwingEq is not yet supported for operation point search)"))
    end

    if ic_guess === nothing
        system_size = systemsize(pg)
        ic_guess = ones(system_size)
    end

    rr = RootRhs_ic(rhs(pg))
    nl_res = nlsolve(rr, ic_guess)
    if converged(nl_res) == true
        return State(pg, nl_res.zero)
    else
        throw(OperationPointError("Failed to find initial conditions on the constraint manifold!"))
    end
end

function check_operationpoint(pg, fixed_point::State; p=nothing, t0=0.)
    ode = rhs(pg)
    M = ode.mass_matrix
    f!(dx, x) = ode(dx, x, p, t0)
    j(x) = (dx = similar(x); jacobian(f!, dx, x))
    λs = eigvals(Array(M), j(fixed_point.vec)) .|> real
    λ = λs[isfinite.(λs)] |> extrema
    println("Jacobian spectrum \nmin : ", first(λ), "\nmax : ",last(λ), "\nstable : ", last(λ)==0)
    return λ
end
check_operationpoint(op::State) = check_operationpoint(op.grid, op.vec)


## ADN global control

function ADN(pg::PowerGrid, DGUnit_Type, P_ref, Q_ref)
    slack_idx = findfirst([node isa SlackAlgebraic for node in pg.nodes])

    Y = PiModel(first(pg.lines))
    sl_func = u -> conj(- Y[1,1] - Y[1, 2]*complex(u[3:4]...))

    args = (pg, slack_idx, P_ref, Q_ref, sl_func)

    ControlledPowerGrid(
        ADN_control,
        pg,
        args...
    )
end

function ADN_control(u, p, t, pg, slack_idx, P_ref, Q_ref, sl_func)

    # this works but raises annoying ND.jl warning every time
    #obs = State(pg, u) # TODO probably better to use graphdata from ND?
    #S_slack = obs[slack_idx, :s]

    # take the high-voltage side as the measurement reference
    # S will be positive since node currents are positive by convention

    S_slack = sl_func(u)

    ΔP = P_ref(t) - real(S_slack)
    ΔQ = Q_ref(t) - imag(S_slack)

    # return ND parameter tuple
    # if p isa Nothing
    #     p_cont = ((ΔP, ΔQ), nothing)
    # else
    #     node_pars, edge_pars = p
    #
    #     p_cont = (
    #         [(ΔP, ΔQ, par...) for par in node_pars],
    #         edge_pars,
    #     )
    # end
    return ((ΔP, ΔQ), nothing)
end
