import Base: @__doc__
using PowerDynamics: AbstractLine
import PowerDynamics: construct_vertex, dimension, symbolsof, showdefinition, AbstractNode
using LinearAlgebra: Diagonal, UniformScaling
using NetworkDynamics: ODEVertex

##

function total_current_real(e_s, e_d)
    # Keeping with the convention of negative sign for outging current
    real_current = 0.0
    imag_current = 0.0
    for e in e_s
        real_current -= e[1]
        imag_current -= e[2]
    end
    for e in e_d
        real_current += e[3]
        imag_current += e[4]
    end
    [real_current; imag_current]
end

##
begin
    Base.@__doc__ struct SlackReal <: AbstractNode
            U
            Y_n
        end
    SlackReal(; U, Y_n = 0) = SlackReal(U, Y_n)
    function construct_vertex(par::SlackReal)
        Ur = real(par.U)
        Ui = imag(par.U)
        Y_n = par.Y_n
        function rhs!(dx, x, e_s, e_d, p, t)
            try
                dx[1] = x[1] - Ur
                dx[2] = x[2] - Ui
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(f! = rhs!, dim = 2, mass_matrix = [0, 0], sym = [:u_r, :u_i])
    end
    symbolsof(::SlackReal) = begin
            [:u_r, :u_i]
        end
    dimension(::SlackReal) = begin
            2
        end
end

##

import PowerDynamics: construct_edge

begin
    @__doc__ struct LineReal <: AbstractLine
            from
            to
            Y
            LineReal(; from, to, Y) = new(from, to, Y)
        end
    function construct_edge(par::LineReal)
        from = par.from
        to = par.to
        G = real(par.Y)
        B = imag(par.Y)
        function rhs!(e, v_s, v_d, p, t)
            # source_voltage = v_s[1] + v_s[2] * im
            # destination_voltage = v_d[1] + v_d[2] * im
            # complex_current = Y * (destination_voltage - source_voltage)
            # current_vector = [complex_current, complex_current]
            e[1] = G * (v_d[1] - v_s[1]) - B * (v_d[2] - v_s[2])#real
            e[2] = G * (v_d[2] - v_s[2]) + B * (v_d[1] - v_s[1])#imag
            e[3] = G * (v_d[1] - v_s[1]) - B * (v_d[2] - v_s[2])#real
            e[4] = G * (v_d[2] - v_s[2]) + B * (v_d[1] - v_s[1])#imag
        end
        return StaticEdge(f! = rhs!, dim = 4)
    end
end

##

begin
    Base.@__doc__ struct PQReal <: AbstractNode
            P
            Q
            Y_n
        end
    PQReal(; P, Q, Y_n = 0) = PQReal(P, Q, Y_n)
    function construct_vertex(par::PQReal)
        P = par.P
        Q = par.Q
        G_n = real(par.Y_n)
        B_n = imag(par.Y_n)
        function rhs!(dx, x, e_s, e_d, p, t)
            i = total_current_real(e_s, e_d) 
            i[1] += G_n * x[1] 
            i[2] += B_n * x[2] 
            p = x[1] * i[1] + x[2] * i[2] 
            q = x[2] * i[1] - x[1] * i[2] 
            #s = u * conj(i)
            #du = (P + (1im) * Q) - s
            try
                dx[1] = P - p #real(du)
                dx[2] = Q - q #imag(du)
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(f! = rhs!, dim = 2, mass_matrix = [0, 0], sym = [:u_r, :u_i])
    end
    symbolsof(::PQReal) = begin
            [:u_r, :u_i]
        end
    dimension(::PQReal) = begin
            2
        end
end

##


begin
    Base.@__doc__ struct ThirdOrderReal <: AbstractNode
            H
            P
            D
            Ω
            E_f
            T_d_dash
            X_q_dash
            X_d_dash
            X_d
            Y_n
        end
    ThirdOrderReal(; H, P, D, Ω, E_f, T_d_dash, X_q_dash, X_d_dash, X_d, Y_n = 0) = ThirdOrderReal(H, P, D, Ω, E_f, T_d_dash, X_q_dash, X_d_dash, X_d, Y_n)
    function construct_vertex(par::ThirdOrderReal)
        H = par.H
        P = par.P
        D = par.D
        Ω = par.Ω
        E_f = par.E_f
        T_d_dash = par.T_d_dash
        X_q_dash = par.X_q_dash
        X_d_dash = par.X_d_dash
        X_d = par.X_d
        G_n = real(par.Y_n)
        B_n = imag(par.Y_n)
        @assert H > 0 "inertia (H) should be >0"
        @assert D >= 0 "damping (D) should be >=0"
        @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
        @assert X_q_dash >= 0
        @assert X_d_dash >= 0 "transient reactance of d-axis (X_d_dash) should be >=0"
        @assert X_d >= 0 "reactance of d-axis (X_d_dash) should be >=0"
        Ω_H = Ω / (2H)
        function rhs!(dx, x, e_s, e_d, p, t)
            i = total_current_real(e_s, e_d) 
            i[1] += G_n * x[1] 
            i[2] += B_n * x[2]
            θ = x[3]
            ω = x[4]
            #i_c = (1im) * i * exp((-1im) * θ)
            #e_c = (1im) * u * exp((-1im) * θ)
            #p = real(u * conj(i))
            p = x[1] * i[1] + x[2] * i[2] 
            e_d = 0
            e_q = x[2] * sin(θ) + x[1] * cos(θ)
            i_d = i[1] * sin(θ) - i[2] * cos(θ)
            i_q = i[2] * sin(θ) + i[1] * cos(θ)
            dθ = ω
            dω = (((P - D * ω) - p) - (X_q_dash - X_d_dash) * i_d * i_q) * Ω_H
            de_q = (1 / T_d_dash) * ((E_f - e_q) + i_d * (X_d - X_d_dash))
            #du = de_q * exp((1im) * θ) + u * (1im) * ω
            try
                dx[1] = de_q * cos(θ) - x[2] * ω #real(du)
                dx[2] = de_q * sin(θ) + x[1] * ω #imag(du)
                dx[3] = dθ
                dx[4] = dω
                return nothing
            catch e
                if typeof(e) === UndefVarError
                    throw(NodeDynamicsError("you need to provide $(e.var)"))
                else
                    throw(e)
                end
            end
        end
        ODEVertex(f! = rhs!, dim = 4, mass_matrix = UniformScaling{Bool}(true), sym = [:u_r, :u_i, :θ, :ω])
    end
    symbolsof(::ThirdOrderReal) = begin
            [:u_r, :u_i, :θ, :ω]
        end
    dimension(::ThirdOrderReal) = begin
            4
        end
end