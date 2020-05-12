import Base: @__doc__
using PowerDynamics: @DynamicNode
import PowerDynamics: construct_vertex, dimension, symbolsof, showdefinition, AbstractNode
using LinearAlgebra: Diagonal
using NetworkDynamics: ODEVertex

begin
    @__doc__ struct DGUnitPLL <: AbstractNode
            I_r
            K_pll
            Y_n
        end
    DGUnitPLL(; I_r, K_pll, Y_n) = DGUnitPLL(I_r, K_pll, Y_n)
    function construct_vertex(par::DGUnitPLL)
        I_r = par.I_r
        K_pll = par.K_pll
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            # current from incoming/outgoing edges, nodal shunt
            i = total_current(e_s, e_d) + Y_n * u
            θ = x[3]

            # PLL
            # (x[2] * cos(θ) - x[1] * sin(θ))
            dθ = - K_pll * abs(u) * sin(θ - angle(u))

            du = i - I_r*exp(θ*1im) # = 0 constraint

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dθ
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(f! = rhs!, dim=3, mass_matrix=Diagonal([0, 0, 1]), sym=Symbol[:u_r, :u_i, :θ])
    end
    symbolsof(::DGUnitPLL) = begin
            [:u_r, :u_i, :θ]
        end
    dimension(::DGUnitPLL) = begin
            3
        end
end


begin
    @__doc__ struct DGUnitPLLPQ <: AbstractNode
            I_r
            K_pll
            S_pq
            Y_n
        end
    DGUnitPLLPQ(; I_r, K_pll, S_pq, Y_n) = DGUnitPLLPQ(I_r, K_pll, S_pq, Y_n)
    function construct_vertex(par::DGUnitPLLPQ)
        I_r = par.I_r
        K_pll = par.K_pll
        S_pq = par.S_pq
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            # current from incoming/outgoing edges, nodal shunt, PQ load background
            i = total_current(e_s, e_d) + Y_n * u - conj(S_pq) / conj(u)
            θ = x[3]

            # PLL
            # (x[2] * cos(θ) - x[1] * sin(θ))
            dθ = - K_pll * abs(u) * sin(θ - angle(u))

            du = i - I_r*exp(θ*1im) # = 0 constraint

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dθ
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(f! = rhs!, dim=3, mass_matrix=Diagonal([0, 0, 1]), sym=Symbol[:u_r, :u_i, :θ])
    end
    symbolsof(::DGUnitPLLPQ) = begin
            [:u_r, :u_i, :θ]
        end
    dimension(::DGUnitPLLPQ) = begin
            3
        end
end

begin
    @__doc__ struct DGUnitPLLPQTracking <: AbstractNode
            K_pll
            K_PT1
            T_PT1
            S_pq
            Y_n
        end
    DGUnitPLLPQTracking(; K_pll, K_PT1, T_PT1, S_pq, Y_n) = DGUnitPLLPQTracking(K_pll, K_PT1, T_PT1, S_pq, Y_n)
    function construct_vertex(par::DGUnitPLLPQTracking)
        K_pll = par.K_pll
        K_PT1 = par.K_PT1
        T_PT1 = par.T_PT1
        S_pq = par.S_pq
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            # current from incoming/outgoing edges, nodal shunt, PQ load background
            i = total_current(e_s, e_d) + Y_n * u - conj(S_pq) / conj(u)
            Vamp = abs(u)

            θ = x[3]
            P_err = x[4]
            Q_err = x[5]
            P_g = x[6]
            Q_g = x[7]

            # power flow tracking
            ΔP, ΔQ = p

            dP_err = - ΔP
            dQ_err = - ΔQ
            dP_g = (K_PT1 * P_err / 2. - P_g) / T_PT1
            dQ_g = (K_PT1 * Q_err / 2. - Q_g) / T_PT1

            # current set point adjustment
            I_r = complex(P_g, -Q_g) / Vamp

            # PLL
            # (x[2] * cos(θ) - x[1] * sin(θ))
            dθ = - K_pll * Vamp * sin(θ - angle(u))

            du = i - I_r*exp(θ*1im) # = 0 constraint

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dθ
                dx[4] = dP_err
                dx[5] = dQ_err
                dx[6] = dP_g
                dx[7] = dQ_g
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(f! = rhs!, dim=7, mass_matrix=Diagonal([0, 0, 1, 1, 1, 1, 1]), sym=Symbol[:u_r, :u_i, :θ, :P_err, :Q_err, :P_g, :Q_g])
    end
    symbolsof(::DGUnitPLLPQTracking) = begin
            [:u_r, :u_i, :θ, :P_err, :Q_err, :P_g, :Q_g]
        end
    dimension(::DGUnitPLLPQTracking) = begin
            7
        end
end

using OrdinaryDiffEq: ODEProblem, DiscreteCallback
import PowerDynamics: NodeShortCircuit

function (nsc::NodeShortCircuit)(powergrid::ControlledPowerGrid)
    # Currently this assumes that the lines are PiModels...
    lines = copy(powergrid.lines)
    l_idx = findfirst(l -> (l.from == nsc.node || l.to == nsc.node) && l isa PiModelLine, lines)

    if l_idx == nothing
        @warn "Node needs to be connected to a PiModelLine to implement NodeShortCircuit"
        return nothing
    end

    if lines[l_idx].from == nsc.node
        lines[l_idx] = PiModelLine(;from=lines[l_idx].from, to=lines[l_idx].to, y = lines[l_idx].y, y_shunt_km = 1/nsc.R,  y_shunt_mk = lines[l_idx].y_shunt_mk)
    elseif lines[l_idx].to == nsc.node
        lines[l_idx] = PiModelLine(;from=lines[l_idx].from, to=lines[l_idx].to, y = lines[l_idx].y, y_shunt_km = lines[l_idx].y_shunt_km,  y_shunt_mk = 1/nsc.R)
    end
    ControlledPowerGrid(powergrid.graph, powergrid.nodes, lines, powergrid.controller)
end

function simulate_nsc(nsc::NodeShortCircuit, powergrid::ControlledPowerGrid, x1, timespan)
    @assert first(timespan) <= nsc.sc_timespan[1] "fault cannot begin in the past"
    @assert nsc.sc_timespan[2] <= last(timespan) "fault cannot end in the future"
    nsc_powergrid = nsc(powergrid)

    problem = ODEProblem{true}(rhs(powergrid), x1, timespan)

    function errorState(integrator)
        sol1 = integrator.sol
        x2 = find_valid_initial_condition(nsc_powergrid, sol1[end]) # Jump the state to be valid for the new system.
        println(x2)
        integrator.f = rhs(nsc_powergrid)
        integrator.u = x2
    end

    function regularState(integrator)
        sol2 = integrator.sol
        x3 = find_valid_initial_condition(nsc_powergrid, sol2[end]) # Jump the state to be valid for the new system.
        println(x3)
        integrator.f = rhs(powergrid)
        integrator.u = x3
    end

    cb1 = DiscreteCallback(((u,t,integrator) -> t in nsc.sc_timespan[1]), errorState)
    cb2 = DiscreteCallback(((u,t,integrator) -> t in nsc.sc_timespan[2]), regularState)
    sol = solve(problem, Rodas4(), callback = CallbackSet(cb1, cb2), tstops=[nsc.sc_timespan[1];nsc.sc_timespan[2]])
    return sol
end
