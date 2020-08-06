using PowerDynamics
using ForwardDiff
using Distributions
import PowerDynamics:rhs
using DifferentialEquations
using LinearAlgebra
using Optim

include("../Random-Pertubation-on-Constrained-Manifolds/ModelingToolKitWrapper.jl")
include("../PDpatches.jl")
include("system.jl")
include("PowerFlow.jl")

"""
# neuer text?
A method to find random pertubations in the case of systems with algebraic
constraints. There is then a set of valid initial conditions which lies on a so
called constrained manifold.
The system must still lie on this manifold after the random pertubation.
This could be applied to Basin Stability Methods in the case of complexer
Generator Methods which contain algebraic constraints.
"""
function RandPertWithConstrains(g::Function, symbol_list,z, op::State, node::Int,
                                interval_V, interval_θ,
                                dist::UnionAll; ProjectionMethod = "Optim",
                                nsteps = 1)

    Jacg,vars = jacobian_wrapper(g, symbol_list) # Getting the symbolic Jacobian
    for steps in 1:nsteps
        Frand = [rand(dist(interval_θ[1],interval_θ[2])), rand(dist(interval_V[1],interval_V[2]))]

        Fstep = Frand / nsteps

        #Jacg = ForwardDiff.gradient(x -> real(g(x[1],x[2],x[3])), z) .+ 1im .* ForwardDiff.gradient(x -> imag(g(x[1],x[2],x[3])), z)

        Jacg_eval = transpose(insert_value(Jacg, vars, z))

        kerJacg = nullspace(Matrix(Jacg_eval')) # kernel of the Jacobian, returns orthogonal vectors

        dz = kerJacg * Fstep

        if ProjectionMethod == "Optim"
            z = ProjOptim(dz, g) # Calculates the projection onto M using Optim

        elseif ProjectionMethod == "ODE"
            z = ProjODE(g, node, op, dz, Q)

        else
            ArgumentError("Please use either Optim or ODE as the ProjectionMethod.")
        end
    end

    return z
end

"""
Calculates and returns the orthogonal Projection of a vector w on the vector v.
"""
function Projvw(v,w)
    dot(w,v) * v * (1 / norm(v)^2)
end

"""
Calculates the Projection from the Tangetial Room TzM to the the Manifold M by
minimizing the function:
    f(x) = a * |dz-p''| + b * |g(p')|.
Where g(p') is the Constraint at point p'.
This function looks so weird because Optim does not support Complex Numbers so
far. See here:
https://stackoverflow.com/questions/61770030/comlex-valued-parameter-estimation-in-julia-using-differentialequations

Inputs:
        dz: The calculated pertubation inside of TzM
        g: The constraint function
        Q: Reactive Power at the perturbed node
Outputs:
        p_dash: Solution of the ODE problem, a Pertubation on M
"""
function ProjOptim(dz::Array{Complex64,1}, g::Function; a = 1.0, b = 100.0)
    len_dz = length(dz)

    Proj(x) = b * abs(g((x[1:2:end] +  1im * x[2:2:end])...)) +
              a * norm(dz - (x[1:2:end] +  1im * x[2:2:end]))

    x0 = Vector{Float64}(undef, 2 * len_dz)

    x0[1:2:end] = real(dz[:])
    x0[2:2:end] = imag(dz[:])

    g((x0[1:2:end] + 1im * x0[2:2:end])...)

    proj = optimize(Proj, x0)

    minimum = proj.minimizer

    p_dash = Vector{Float64}(undef, len_dz)

    p_dash[:] = abs.(minimum[1:2:end] + 1im * minimum[2:2:end])

    return p_dash
end

function ProjOptim(dz::Array{Float64,1}, g::Function; a = 1.0, b = 100.0)
    Proj(x) = b * abs(g(x...)) +
              a * norm(dz - x)

    proj = optimize(Proj, dz)

    minimum = proj.minimizer

    return minimum
end

"""
Calculates the Projection from the Tangetial Room TzM to the the Manifold M by
solving a ODE of the form:
    f(x) = y' = - g(dz, y).
Where dz is the result of the pertubation and a point in TzM.
Which could result in a fixed point g(x_v, y*) = 0.
Inputs:
        dz: The calculated pertubation inside of TzM
        g: The constraint function
        Q: Reactive Power at the perturbed node
        endtime: Integration time of the problem
Outputs:
        p_dash: The result of optimization, a Pertubation on M
"""
function ProjODE(g::Function,node,op,dz,Q;endtime = 10.0)
    tspan = (0.0,endtime)

    u0 = [dz[1], dz[2], dz[3], g(dz[1], dz[2], dz[3]), op[node,:ω]]
    vars = powergrid.nodes[node]
    #function constraint(du,u,p,t)
    #    du[4] = - (Q - 3 * u[3] * u[2] * exp(-2im * u[1]))
    #end
    function constraint(du,u,p,t)
        Ω_H = vars.Ω / (2 * vars.H)
        i_c = 1im * u[2] * exp(-1im * u[1])
        e_c = 1im * u[3] * exp(-1im * u[1])
        p = real(u[3] * conj(u[2]))
        e_d = 0                         # simplification from K.Schmietendorf paper
        e_q = imag(e_c)
        i_d = real(i_c)
        i_q = imag(i_c)

        du[1] = u[5]
        du[5] = (vars.P - vars.D * u[5] - p - (vars.X_q_dash - vars.X_d_dash) * i_d * i_q) * Ω_H
        de_q = (1 / vars.T_d_dash) * (vars.E_f - e_q + i_d * (vars.X_d - vars.X_d_dash))
        # -> u = e_q * exp(1im * θ)
        du[3] = de_q * exp(1im * u[1]) + u[3] * 1im * u[5]
        du[4] = u[4] - (Q - 3 * u[3] * u[2] * exp(-2im * u[1]))
    end

    prob = ODEProblem(constraint,u0,tspan)

    sol = solve(prob, Rosenbrock23(autodiff=false), abstol=1e-8,reltol=1e-6,tspan=Inf, force_dtmin = true)
    p_dash = sol[1:3,end]
end

function rhs(pg::PowerGrid)
    network_dynamics(map(construct_vertex, pg.nodes), map(construct_edge, pg.lines), pg.graph)
end

"""
A method to assures that the conditions given from g stay on the contrained
manifold. A DifferentialEquations Callback is used.
This DOES NOT work because of a type missmatch regarding complex numbers
Inputs:
    pg: Power grid, a graph containg nodes and lines
    endtime: Endtime of the simulation
    u0: current state of the powergrid pg
Outputs:
    solution: a PowerGrid solution containg the timeseries and the Powergrid pg
"""
function SimManifoldCd(pg::PowerGrid, Q, u0::State, endtime)
    ode = ODEProblem(rhs(pg), u0.vec, (0.0, endtime), nothing)

    function g(resid, x, e_s, e_d, p, t)
            i = total_current(e_s, e_d) + Y_n * (x[1] + x[2] * im)
            u = x[1] + x[2] * im
            θ = x[3]
            ω = x[4]
            resid[1] = abs(Q) - 3 * abs(i * u * exp((-2im) * θ))
    end

    function g1(resid, x, e_s, e_d, p, t)
        resid[1] = x[2]^2 + x[1]^2 - 2
        resid[2] = 0
    end


    cb = ManifoldProjection(g1)


    dqsol = solve(ode, Rodas5(autodiff = false), save_everystep=false, callback=cb)
    solution = PowerGridSolution(dqsol, pg)
    return solution
end
