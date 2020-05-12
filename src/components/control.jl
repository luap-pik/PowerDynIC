function ADN(pg::PowerGrid, DGUnit_Type, P_ref, Q_ref)
    slack_idx = findfirst([node isa SlackAlgebraic for node in pg.nodes])

    args = (pg, slack_idx, P_ref, Q_ref)

    ControlledPowerGrid(
        ADN_control,
        pg,
        args...
    )
end

function ADN_control(u, p, t, pg, slack_idx, P_ref, Q_ref)

    obs = State(pg, u) # TODO probably better to use graphdata from ND?

    # take the high-voltage side as the measurement reference
    # S will be positive since node currents are positive by convention
    S_slack = obs[slack_idx, :s]

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
