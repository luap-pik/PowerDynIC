function CIGRE_dynamic(DG_unit_type, busses_static, lines, T, elist, Zs, Yshs, k_FRT_b)

    busses_dynamic = copy(busses_static)

    DG_locs = 2:12

    for i in DG_locs
        P_bkgrnd = 0.
        Q_bkgrnd = 0.
        try
            P_bkgrnd = real(busses_static[i].S)
            Q_bkgrnd = imag(busses_static[i].S)
        catch
            P_bkgrnd = busses_static[i].P # always < 0
            Q_bkgrnd = busses_static[i].Q
        end
        busses_dynamic[i] = eval(Symbol(DG_unit_type))(
            P_bkgrnd = P_bkgrnd,
            Q_bkgrnd = Q_bkgrnd,
            Pr_l_min = 0., # MW, Holm: 10e12 / base_power (= not used)
            Pr_l_max = 1.0, #1.0, # MW, Holm: 10e12 / base_power (= not used)
            Qr_l_min = -1.0, # MW, Holm: 10e12 / base_power (= not used)
            Qr_l_max = 1.0, # MW, Holm: 10e12 / base_power (= not used)
            Pr_del_min = 0., #-1.0, # Holm: 0 MW - here -2.0 because overall signs are not clear
            Pr_del_max = 1.0, # Holm: 200 MW (= not used)
            Qr_del_min = -1.0, # Holm: -60 MW (= not used)
            Qr_del_max = 1.0, # Holm: 60 MW (= not used)
            T_PT1 = 0.1, # seconds
            K_PT1 = 1.0, # PT1 gain relevant for fix point
            Vac_ref = lines[1].Uus, #* sqrt(2/3),
            V_dead = 0.1 * lines[1].Uus, #* sqrt(2/3),
            k_FRT = k_FRT_b, #/ base_current,
            imax = 1e6/20e3 / base_current / sqrt(3), # MW/kV
            K_pll = 0.1 * base_voltage,
            Y_shunt = t -> 0.#-0.1 * SmoothStep(t, 3.15, 3.; order=100) #0.,
        )
    end

    begin
        lines_dynamic = Array{AbstractLine,1}([])
        push!(lines_dynamic, T)
        for (e, Z, Ysh) in zip(elist, Zs, Yshs)
            push!(
                lines_dynamic,
                RLLine(
                    from = first(e),
                    to = last(e),
                    R = real(Z) * base_admittance,
                    L = (imag(Z) / ω)  * base_admittance,
                    ω0 = ω
                ),
            )
        end
    end

busses_dynamic, lines_dynamic
end
