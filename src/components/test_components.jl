import Base: @__doc__
using PowerDynamics: @DynamicNode
import PowerDynamics: construct_vertex, dimension, symbolsof, showdefinition, AbstractNode
using LinearAlgebra: Diagonal

begin
    @__doc__ struct CSImin_limit <: AbstractNode
            I_r
            imax
            Y_n
        end
    CSImin_limit(; I_r, imax, Y_n=0) = CSImin_limit(I_r, imax, Y_n)
    function construct_vertex(par::CSImin_limit)
        I_r = par.I_r
        imax = par.imax
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            i = total_current(e_s, e_d) + Y_n * u

            id_r = SmoothClamp(real(i), imax; scale=500.)
            iq_r = SmoothClamp(imag(i), sqrt(abs(imax^2 - id_r^2)); scale=500.)

            du = complex(id_r, iq_r) - I_r

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(f! = rhs!, dim=2, mass_matrix=[0.0, 0.0], sym=Symbol[:u_r, :u_i])
    end
    symbolsof(::CSImin_limit) = begin
            [:u_r, :u_i]
        end
    dimension(::CSImin_limit) = begin
            2
        end
end

begin
    @__doc__ struct CSImin_FRT <: AbstractNode
            I_r
            Vac_ref
            V_dead
            imax
            k_FRT
            Y_n
        end
    CSImin_FRT(; I_r, Vac_ref, V_dead, imax, k_FRT, Y_n=0) = CSImin_FRT(I_r, Vac_ref, V_dead, imax, k_FRT, Y_n)
    function construct_vertex(par::CSImin_FRT)
        I_r = par.I_r
        Vac_ref = par.Vac_ref
        V_dead = par.V_dead
        imax = par.imax
        k_FRT = par.k_FRT
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            i = total_current(e_s, e_d) + Y_n * u

            Va_amp = abs(u)

            #err_state = SmoothStep(Va_amp, Vac_ref + V_dead, Vac_ref - V_dead; order=100.)
            err_state = Vac_ref - V_dead < Va_amp < Vac_ref + V_dead

            #fault err != 0
            id_r = SmoothClamp(real(i), imax; scale=200.)
            iq_r = SmoothClamp(imag(i), sqrt(abs(imax^2 - id_r^2)); scale=200.)

            iq_plus = 0.
            if Va_amp > Vac_ref + V_dead
                iq_plus = - k_FRT * (Va_amp - (Vac_ref + V_dead))
            elseif Va_amp < Vac_ref - V_dead
                iq_plus = - k_FRT * (Va_amp - (Vac_ref - V_dead))
            end
            iq_r_plus = SmoothClamp(iq_plus, imax; scale=200.)
            id_r_plus = SmoothClamp(iq_r_plus, sqrt(abs(imax^2 - iq_r^2)); scale=200.)

            #i_mix = err_state * complex(id_r, iq_r) + (1. - err_state) * complex(id_r_plus, iq_r_plus)
            i_mix = err_state ? complex(id_r, iq_r) : complex(id_r_plus, iq_r_plus)

            du = i_mix - I_r

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(f! = rhs!, dim=2, mass_matrix=[0.0, 0.0], sym=Symbol[:u_r, :u_i])
    end
    symbolsof(::CSImin_FRT) = begin
            [:u_r, :u_i]
        end
    dimension(::CSImin_FRT) = begin
            2
        end
end

begin
    @__doc__ struct CSImin_PLL <: AbstractNode
            I_r
            K_pll
            Y_n
        end
    CSImin_PLL(; I_r, K_pll, Y_n=0) = CSImin_PLL(I_r, K_pll, Y_n)
    function construct_vertex(par::CSImin_PLL)
        I_r = par.I_r
        K_pll = par.K_pll
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            u = complex(x[1], x[2])
            i = total_current(e_s, e_d) + Y_n * u
            theta = x[3]

            dtheta = - K_pll * abs(u) * sin(theta - angle(u))
            du = i - I_r

            try
                dx[1] = real(du)
                dx[2] = imag(du)
                dx[3] = dtheta
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
    symbolsof(::CSImin_PLL) = begin
            [:u_r, :u_i, :θ]
        end
    dimension(::CSImin_PLL) = begin
            3
        end
end
