using PowerDynamics: AbstractLine, construct_edge
using LinearAlgebra: I
import PowerDynamics: PiModel

begin
    @__doc__ struct RLLine <: AbstractLine
            from
            to
            R
            L
            ω0
            RLLine(; from, to, R, L, ω0) = new(from, to, R, L, ω0)
        end
    function PowerDynamics.construct_edge(par::RLLine)
        from = par.from
        to = par.to
        R = par.R
        L = par.L
        ω0 = par.ω0
        @assert R > 0 "resistance (R) should be positive"
        @assert L > 0 "inductance (L) should be positive"
        @assert ω0 > 0 "rated frequency (ω0) should be positive"
        Z = [R -ω0*L; ω0*L R]
        function rhs!(de, e, v_s, v_d, p, t)
            # the different minus signs are due to the PowerDynamics sign convention for currents
            i_left_right = [e[1]; e[2]]
            i_right_left = [e[3]; e[4]]
            v_left = [v_s[1]; v_s[2]]
            v_right = [v_d[1]; v_d[2]]
            v_left_right = (- Z * i_left_right .- v_left .+ v_right) ./ L
            v_right_left = (- Z * i_right_left .+ v_right .- v_left) ./ L
            de .= [v_left_right; v_right_left]
        end
        #TODO check if there's a bug in the macro missing the empty space
        return ODEEdge(f! = rhs!, dim=4, mass_matrix=I, sym=Symbol[:i_r, :i_i, :ir_r, :ir_i])
    end
    symbolsof(::RLLine) = begin
            [:i_r, :i_i, :ir_r, :ir_i]
    end
    dimension(::RLLine) = begin
            4
    end
end

PiModel(l::RLLine) = PiModel(1/(l.R + im *l.ω0*l.L), 0, 0, 1, 1)
