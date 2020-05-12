function calcY(line)
    t_km = t_mk = 1
    y = line.y
    y_shunt_km = line.y_shunt_km
    y_shunt_mk = line.y_shunt_mk
    Π = zeros(Complex{Float64}, 2, 2)
    Π[1, 1] = abs2(t_km) * (y + y_shunt_km) # Our sign convention is opposite for the source of the edge
    Π[1, 2] = - conj(t_km) * t_mk * y # Our sign convention is opposite for the source of the edge
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
end

function map_complex_2D(array)
    n = length(array)
    out = zeros(2n)
    for i = 1:n
        out[2i-1] = real(array[i])
        out[2i] = imag(array[i])
    end
    return out
end

function map_2D_complex(array)
    return array[1:2:end] .+ 1im .* array[2:2:end]
end

function niceprint(pf; digits=2)
    v_angles = v_angle(pf; digits=2)
    for val in v_angles
        println(val)
    end
end

function v_angle(pf; digits=2)
    voltage_solution = abs.(map_2D_complex(pf)) .|> x -> round(x; digits=digits)
    angle_solution = angle.(map_2D_complex(pf)) .|> rad2deg .|> x -> round(x; digits=digits)
    return [string(v)*"∠"*string(a) for (v, a) in zip(voltage_solution, angle_solution)]
end


struct rr
    #=
    mpm der Projektionsoperator auf das Bild der MAssmatrix. in diesem bild
    mus f(x) liegen damit Mx' = f(x) eine lösung hat.
    =#
    rhs
    mpm
    p
end

function (rr::rr)(f_x, x)
    rr.rhs(f_x, x, rr.p, 0.0)
    f_x .-= rr.mpm * f_x
end

#TODO: fix
function determine_valid_ic(rhs, p, guess)
    # pfsol = pf_sol(busses_static, lines; matfile="$dir/LFCigre.mat", outfile = nothing)
    # idx = [1, 2]
    # for v_idx in cpg.volt_idxs
    #     ic_guess[v_idx] .= pf[idx] # in kV
    #     idx .+= 2
    # end
    # #set the PLL to the voltage angle solution
    # θ = deg2rad.(power_flow_sol.A)
    # θ_idx = [first(d_idx) for d_idx in cpg.dg_unit_idxs]
    # for (i, idx) in enumerate(θ_idx)
    #     ic_guess[first(idx)] = θ[i]
    # end

    mpm = Diagonal(rhs.mass_matrix) * pinv(Diagonal(rhs.mass_matrix))
    pg_rr = rr(rhs, mpm, p)
    nlsol = nlsolve(pg_rr, guess; xtol=1E-12) # autodiff = :forward
    println("nlsolve converged ", nlsol.f_converged)

    dx = similar(nlsol.zero)
    pg_rr = rr(rhs, mpm, p)
    pg_rr(dx, nlsol.zero)
    println(sum(abs, dx))

    return nlsol.f_converged ? nlsol.zero : guess
end

#TODO: call this function until returns true
#TODO: PiModel as function of t
function determine_random_valid_ic(rhs, p, guess, ϵ)
    mpm = Diagonal(rhs.mass_matrix) * pinv(Diagonal(rhs.mass_matrix))
    pg_rr = rr(rhs, mpm, p)
    nlsol = nlsolve(pg_rr, guess .+ ϵ .* randn(length(guess)); autodiff = :forward, xtol=1E-2)
    println("nlsolve converged ", nlsol.f_converged)

    dx = similar(nlsol.zero)
    pg_rr = rr(rhs, mpm, p)
    pg_rr(dx, nlsol.zero)
    println(sum(abs, dx))

    return nlsol.f_converged, nlsol.f_converged ? nlsol.zero : guess
end
