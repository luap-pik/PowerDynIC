import Base: @__doc__
using PowerDynamics: @DynamicNode
import PowerDynamics:
    construct_vertex, dimension, symbolsof, showdefinition, AbstractNode
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
            dθ = -K_pll * abs(u) * sin(θ - angle(u))

            du = i - I_r * exp(θ * 1im) # = 0 constraint

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
        ODEVertex(
            f! = rhs!,
            dim = 3,
            mass_matrix = Diagonal([0, 0, 1]),
            sym = Symbol[:u_r, :u_i, :θ],
        )
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
            dθ = -K_pll * abs(u) * sin(θ - angle(u))

            du = i - I_r * exp(θ * 1im) # = 0 constraint

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
        ODEVertex(
            f! = rhs!,
            dim = 3,
            mass_matrix = Diagonal([0, 0, 1]),
            sym = Symbol[:u_r, :u_i, :θ],
        )
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
    DGUnitPLLPQTracking(; K_pll, K_PT1, T_PT1, S_pq, Y_n) =
        DGUnitPLLPQTracking(K_pll, K_PT1, T_PT1, S_pq, Y_n)
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
            P_int = x[4]
            Q_int = x[5]
            P_g = x[6]
            Q_g = x[7]

            # power flow tracking
            ΔP, ΔQ = p

            dP_int = -ΔP
            dQ_int = -ΔQ
            dP_g = (K_PT1 * P_int / 2.0 - P_g) / T_PT1
            dQ_g = (K_PT1 * Q_int / 2.0 - Q_g) / T_PT1

            # current set point adjustment
            I_r = complex(P_g, -Q_g) / Vamp

            # PLL
            # (x[2] * cos(θ) - x[1] * sin(θ))
            dθ = -K_pll * Vamp * sin(θ - angle(u))

            du = i - I_r * exp(θ * 1im) # = 0 constraint

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dθ
                dx[4] = dP_int
                dx[5] = dQ_int
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
        ODEVertex(
            f! = rhs!,
            dim = 7,
            mass_matrix = Diagonal([0, 0, 1, 1, 1, 1, 1]),
            sym = Symbol[:u_r, :u_i, :θ, :P_int, :Q_int, :P_g, :Q_g],
        )
    end
    symbolsof(::DGUnitPLLPQTracking) = begin
        [:u_r, :u_i, :θ, :P_int, :Q_int, :P_g, :Q_g]
    end
    dimension(::DGUnitPLLPQTracking) = begin
        7
    end
end

begin
    @__doc__ struct DGUnitTrackingFRT <: AbstractNode
        K_pll
        K_PT1
        T_PT1
        K_FRT
        I_max
        Vref
        Vdead
        S_pq
        Y_n
    end
    DGUnitTrackingFRT(;
        K_pll,
        K_PT1,
        T_PT1,
        K_FRT,
        I_max,
        Vref,
        Vdead,
        S_pq,
        Y_n,
    ) = DGUnitTrackingFRT(
        K_pll,
        K_PT1,
        T_PT1,
        K_FRT,
        I_max,
        Vref,
        Vdead,
        S_pq,
        Y_n,
    )
    function construct_vertex(par::DGUnitTrackingFRT)
        K_pll = par.K_pll
        K_PT1 = par.K_PT1
        T_PT1 = par.T_PT1
        K_FRT = par.K_FRT
        I_max = par.I_max
        Vref = par.Vref
        Vdead = par.Vdead
        S_pq = par.S_pq
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            # current from incoming/outgoing edges, nodal shunt, PQ load background
            i = total_current(e_s, e_d) + Y_n * u - conj(S_pq) / conj(u)
            Vamp = abs(u)

            θ = x[3]
            P_int = x[4]
            Q_int = x[5]
            P_g = x[6]
            Q_g = x[7]

            # power flow tracking
            ΔP, ΔQ = p

            # detect error state
            over_voltage = Vamp > Vref + Vdead
            under_voltage = Vamp > Vref - Vdead
            smooth_fault_state(x) =
                SmoothStep(x, Vref + Vdead, Vref - Vdead; order = 100)

            # don't integrate the error during faults
            dP_int = -ΔP * smooth_fault_state(Vamp)
            dQ_int = -ΔQ * smooth_fault_state(Vamp)

            dP_g = (K_PT1 * P_int / 2.0 - P_g) / T_PT1
            dQ_g = (K_PT1 * Q_int / 2.0 - Q_g) / T_PT1

            # current set point adjustment
            I_r = complex(P_g, -Q_g) / Vamp

            # FRT
            # additional reactive current during fault
            iqplus = 0.0
            if over_voltage
                println("over")
                iqplus = -K_FRT * (Vamp - (Vref + Vdead))
            elseif under_voltage
                println("under")
                iqplus = -K_FRT * (Vamp - (Vref - Vdead))
            end

            I_r += complex(0.0, -iqplus)


            # PLL
            # (x[2] * cos(θ) - x[1] * sin(θ))
            dθ = -K_pll * Vamp * sin(θ - angle(u))

            du = i - I_r * exp(θ * 1im) # = 0 constraint

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dθ
                dx[4] = dP_int
                dx[5] = dQ_int
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
        ODEVertex(
            f! = rhs!,
            dim = 7,
            mass_matrix = Diagonal([0, 0, 1, 1, 1, 1, 1]),
            sym = Symbol[:u_r, :u_i, :θ, :P_int, :Q_int, :P_g, :Q_g],
        )
    end
    symbolsof(::DGUnitTrackingFRT) = begin
        [:u_r, :u_i, :θ, :P_int, :Q_int, :P_g, :Q_g]
    end
    dimension(::DGUnitTrackingFRT) = begin
        7
    end
end

begin
    @__doc__ struct StaticGenerator <: AbstractNode
        K_pll
        S_ref
        K_FRT
        I_max
        Vref
        Vdead
        Y_n
    end
    StaticGenerator(;
        K_pll,
        S_ref,
        K_FRT,
        I_max,
        Vref,
        Vdead,
        Y_n,
    ) = StaticGenerator(
        K_pll,
        S_ref,
        K_FRT,
        I_max,
        Vref,
        Vdead,
        Y_n,
    )
    function construct_vertex(par::StaticGenerator)
        K_pll = par.K_pll
        S_ref = par.S_ref
        K_FRT = par.K_FRT
        I_max = par.I_max
        Vref = par.Vref
        Vdead = par.Vdead
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            # current from incoming/outgoing edges, nodal shunt, PQ load background
            i = total_current(e_s, e_d) + Y_n(t) * u
            Vamp = abs(u)

            θ = x[3]
            i_d_r = x[4]
            i_q_r = x[5]

            # PLL
            dθ = -K_pll * Vamp * sin(θ - angle(u))

            # FRT
            # additional reactive current during fault
            if Vamp > Vref + Vdead
                iqplus = K_FRT * (Vamp - (Vref + Vdead)) # neg for over_voltage
            elseif Vamp < Vref - Vdead
                iqplus = K_FRT * (Vamp - (Vref - Vdead)) # pos for under_voltage
            else
                iqplus = 0.
            end

            # current set point adjustment
            I_r = conj(S_ref) / conj(u) #exp(θ * 1im) *

            # save values for plotting
            di_d_r = i_d_r - real(I_r)
            di_q_r = i_q_r - imag(I_r)

            # current limiter
            if Vamp > Vref + Vdead || Vamp < Vref - Vdead
                iq_fault = Limit(imag(I_r) + iqplus, -I_max, I_max)
                id_lim = sqrt(I_max^2 - iq_fault^2)
                id_fault = Limit(real(I_r), 0., id_lim) #-id_lim, id_lim)
                I_r_limited = complex(id_fault, iq_fault)
            else
                id_normal = Limit(real(I_r), 0., I_max) #-I_max, I_max)
                iq_lim = sqrt(I_max^2 - id_normal^2)
                iq_normal = Limit(imag(I_r), -iq_lim, iq_lim)
                I_r_limited = complex(id_normal, iq_normal)
            end

            du = i - I_r_limited # = 0 constraint

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dθ
                dx[4] = di_d_r
                dx[5] = di_q_r
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(
            f! = rhs!,
            dim = 5,
            mass_matrix = Diagonal([0, 0, 1, 0, 0]),
            sym = Symbol[:u_r, :u_i, :θ, :i_d_r, :i_q_r],
        )
    end
    symbolsof(::StaticGenerator) = begin
        [:u_r, :u_i, :θ, :i_d_r, :i_q_r]
    end
    dimension(::StaticGenerator) = begin
        5
    end
end

function Limit(x, l, u)
    if x < l
        return l
    elseif x > u
        return u
    else
        return x
    end
end



begin
    @__doc__ struct StaticGeneratorObs <: AbstractNode
        K_pll
        S_ref
        K_FRT
        I_max
        Vref
        Vdead
        Y_n
    end
    StaticGeneratorObs(; K_pll, S_ref, K_FRT, I_max, Vref, Vdead, Y_n) =
        StaticGeneratorObs(K_pll, S_ref, K_FRT, I_max, Vref, Vdead, Y_n)
    function construct_vertex(par::StaticGeneratorObs)
        K_pll = par.K_pll
        S_ref = par.S_ref
        K_FRT = par.K_FRT
        I_max = par.I_max
        Vref = par.Vref
        Vdead = par.Vdead
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            # current from incoming/outgoing edges, nodal shunt, PQ load background
            i = total_current(e_s, e_d) + Y_n * u
            Vamp = abs(u)

            # this is the only dynamic variable
            θ = x[3]

            # to keep track of the internal currents for debugging ...
            i_d_r = x[4] # behind FRT
            i_q_r = x[5] # behind FRT
            i_d_r_b = x[6] # before FRT
            i_q_r_b = x[7] # before FRT
            iqp = x[8]

            # PLL
            dθ = -K_pll * Vamp * sin(θ - angle(u))

            # FRT
            # additional reactive current during fault
            if Vamp > Vref + Vdead
                iqplus = K_FRT * (Vamp - (Vref + Vdead)) # pos for over_voltage
            elseif Vamp < Vref - Vdead
                iqplus = K_FRT * (Vamp - (Vref - Vdead)) # neg for under_voltage
            else
                iqplus = 0.0
            end

            # current set point
            I_r = conj(S_ref) / Vamp

            # save values for plotting
            di_d_r_b = i_d_r_b - real(I_r)
            di_q_r_b = i_q_r_b - imag(I_r)
            diqp = iqp - iqplus

            # FRT limiter
            if Vamp > Vref + Vdead || Vamp < Vref - Vdead
                iq_fault = Limit(imag(I_r) + iqplus, -I_max, I_max)
                id_lim = sqrt(I_max^2 - iq_fault^2)
                id_fault = Limit(real(I_r), 0.0, id_lim) #-id_lim, id_lim)
                I_r_limited = complex(id_fault, iq_fault)
            else
                id_normal = Limit(real(I_r), 0.0, I_max) #-I_max, I_max)
                iq_lim = sqrt(I_max^2 - id_normal^2)
                iq_normal = Limit(imag(I_r), -iq_lim, iq_lim)
                I_r_limited = complex(id_normal, iq_normal)
            end

            # save values for plotting
            di_d_r = i_d_r - real(I_r_limited)
            di_q_r = i_q_r - imag(I_r_limited)

            du = i - exp(θ * 1im) * I_r_limited # = 0 constraint

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dθ
                dx[4] = di_d_r
                dx[5] = di_q_r
                dx[6] = di_d_r_b
                dx[7] = di_q_r_b
                dx[8] = diqp
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(
            f! = rhs!,
            dim = 8,
            mass_matrix = Diagonal([0, 0, 1, 0, 0, 0, 0, 0]),
            sym = Symbol[
                :u_r,
                :u_i,
                :θ,
                :i_d_r,
                :i_q_r,
                :i_d_r_b,
                :i_q_r_b,
                :iqp,
            ],
        )
    end
    symbolsof(::StaticGeneratorObs) = begin
        [:u_r, :u_i, :θ, :i_d_r, :i_q_r, :i_d_r_b, :i_q_r_b, :iqp]
    end
    dimension(::StaticGeneratorObs) = begin
        8
    end
end

# x = -2:.01:2
# plot(x, SmoothStep.(x, 1, -1))
# plot!(x, SmoothClamp.(x, 1, -1))
function smooth_step(x, loc; order = 10)
    (1 .+ tanh.(order * (x - loc))) / 2.0
end

function smooth_cutoff(value, threshold; scale=10.)
    log(exp(scale * value) + exp.(scale * threshold)) / scale
end

begin
    @__doc__ struct DGUnit <: AbstractNode
        K_pll
        K_PT1
        T_PT1
        K_FRT
        I_max
        Vref
        Vdead
        S_pq
        Y_n
    end
    DGUnit(; K_pll, K_PT1, T_PT1, K_FRT, I_max, Vref, Vdead, S_pq, Y_n) =
        DGUnit(K_pll, K_PT1, T_PT1, K_FRT, I_max, Vref, Vdead, S_pq, Y_n)
    function construct_vertex(par::DGUnit)
        K_pll = par.K_pll
        K_PT1 = par.K_PT1
        T_PT1 = par.T_PT1
        K_FRT = par.K_FRT
        I_max = par.I_max
        Vref = par.Vref
        Vdead = par.Vdead
        S_pq = par.S_pq
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            # current from incoming/outgoing edges, nodal shunt, PQ load background
            i = total_current(e_s, e_d) + Y_n * u - conj(S_pq) / conj(u)
            Vamp = abs(u)

            # this are the dynamic variables
            θ = x[3]
            P_int = x[4]
            Q_int = x[5]
            P_g = x[6]
            Q_g = x[7]

            # to keep track of the internal currents for debugging ...
            iqp = x[8]

            # power flow tracking
            ΔP, ΔQ = p

            # detect error state
            over_voltage = Vamp > Vref + Vdead
            under_voltage = Vamp < Vref - Vdead
            # smooth_fault_state(x) = SmoothStep(x, Vref + Vdead, Vref - Vdead; order = 100)

            # don't integrate the error during faults
            if over_voltage | under_voltage
                dP_int = 0.
                dQ_int = 0.
                dP_g = 0.
                dQ_g = 0.
            else
                dP_int = -ΔP #* smooth_fault_state(Vamp)
                dQ_int = -ΔQ #* smooth_fault_state(Vamp)
                dP_g = (K_PT1 * P_int / 2.0 - P_g) / T_PT1
                dQ_g = (K_PT1 * Q_int / 2.0 - Q_g) / T_PT1
            end

            # dP_int = -ΔP
            # dQ_int = -ΔQ

            # dP_g = (K_PT1 * P_int / 2.0 - P_g) / T_PT1
            # dQ_g = (K_PT1 * Q_int / 2.0 - Q_g) / T_PT1

            # PLL
            dθ = -K_pll * Vamp * sin(θ - angle(u))

            # FRT
            # additional reactive current during fault
            if over_voltage
                iqplus = K_FRT * (Vamp - (Vref + Vdead)) # pos for over_voltage
            elseif under_voltage
                iqplus =-K_FRT * ((Vref - Vdead) - Vamp) # neg for under_voltage
            else
                iqplus = 0.0
            end

            # current set point
            I_r = conj(complex(P_g, Q_g)) / Vamp

            # save values for plotting
            diqp = iqp - iqplus

            # FRT limiter
            if over_voltage | under_voltage
                iq_fault = Limit(imag(I_r) + iqplus, -I_max, I_max)
                id_lim = sqrt(I_max^2 - iq_fault^2)
                id_fault = Limit(real(I_r), 0.0, id_lim) #-id_lim, id_lim)
                I_r_limited = complex(id_fault, iq_fault)
            else
                id_normal = Limit(real(I_r), 0.0, I_max) #-I_max, I_max)
                iq_lim = sqrt(I_max^2 - id_normal^2)
                iq_normal = Limit(imag(I_r), -iq_lim, iq_lim)
                I_r_limited = complex(id_normal, iq_normal)
            end

            du = i - exp(θ * 1im) * I_r_limited  # = 0 constraint

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dθ
                dx[4] = dP_int
                dx[5] = dQ_int
                dx[6] = dP_g
                dx[7] = dQ_g
                dx[8] = diqp
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(
            f! = rhs!,
            dim = 8,
            mass_matrix = Diagonal([0, 0, 1, 1, 1, 1, 1, 0]),
            sym = Symbol[
                :u_r,
                :u_i,
                :θ,
                :P_int,
                :Q_int,
                :P_g,
                :Q_g,
                :iqp,
            ],
        )
    end
    symbolsof(::DGUnit) = begin
        [:u_r, :u_i, :θ, :P_int, :Q_int, :P_g, :Q_g, :iqp]
    end
    dimension(::DGUnit) = begin
        8
    end
end


function clamp(x, lowerlimit, upperlimit)
    if (x < lowerlimit)
        x = lowerlimit
    elseif (x > upperlimit)
        x = upperlimit
    end
    return x
end
function smootherstep(edge0, edge1, x)
    # Scale, and clamp x to 0..1 range
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.)
    # Evaluate polynomial
    return x * x * x * (x * (x * 6 - 15) + 10)
end
