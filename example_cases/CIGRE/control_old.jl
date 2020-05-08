
mutable struct Controlled_PG
    controller
    open_loop_dyn
    dg_unit_idxs
    volt_idxs
    slack_connect
    P_ref # function of t
    Q_ref # function of t
end
function (cont::Controlled_PG)(du, u, p, t)
    # Get controller state
    # p holds the initial error state of the voltage limit
    p_cont = cont.controller(u, p, t, cont.volt_idxs, cont.slack_connect, cont.P_ref, cont.Q_ref)

    # Calculate the derivatives
    cont.open_loop_dyn(du, u, p_cont, t)
end

function Controlled_PG(pg, DGUnit_Type, P_ref, Q_ref)
    volt_idxs = Array{Array{Int64,1},1}() # real and imaginary part of voltage
    dg_unit_idxs = Array{Array{Int64,1},1}() # internal DGUnit indices

    offset = 1
    for (idx, node) in enumerate(pg.nodes)
        dim = dimension(node)
        push!(volt_idxs, offset:(offset+1))
        if node isa DGUnit_Type
            push!(dg_unit_idxs, (offset+2):(offset+dim-1))
        end
        offset += dim
    end

    slack_connect_idx = zeros(2)
    for line in pg.lines
        if line.from == 1
            slack_connect_idx = volt_idxs[line.to]
        end
    end

    transformer = pg.lines[1]
    Π = PiModel(transformer)

    slack_connect = slack_connect_idx, Π
    Controlled_PG(
        control_loop,
        rhs(pg),
        dg_unit_idxs,
        volt_idxs,
        slack_connect,
        P_ref,
        Q_ref,
    )
end

function control_loop(u, p, t, volt_idxs, slack_connect, P_ref, Q_ref)
    u_s = [complex(u[1], u[2]); complex(u[slack_connect[1][1]], u[slack_connect[1][2]])]
    Π = slack_connect[2]
    S_slack = u_s .* conj(Π * u_s)

    # take the high-voltage side as measurement reference
    # The minus sign is important and due to the sign convention in PowerDynamics
    # within the definition of PiModel. Π ↦ Π .* [-1 -1; 1 1]

    P_meas = real(-S_slack[1])
    Q_meas = imag(-S_slack[1])
    ΔP = P_ref(t) - P_meas # P_ref in Mw
    ΔQ = Q_ref(t) - Q_meas # Q_ref in Mw

    # return ND parameter tuple
    node_pars, edge_pars = p
    ([(ΔP, ΔQ, global_limit, local_limit) for (global_limit, local_limit) in node_pars], nothing)
end

############### callback for local voltage_limit_achieved #############

struct LocalVoltageLimitCondition
    UsollMV
    vlim
    volt_idxs
end
function (vla::LocalVoltageLimitCondition)(out, u, t, integrator) # condition
    ΔV = [norm(u[v_idx]) for v_idx in vla.volt_idxs[2:end]] .- vla.UsollMV
    out .= vla.vlim .- vla.UsollMV .- ΔV # bitte nicht mehr ändern!
    #return nothing
end
function EnterLocalErrorState(integrator, event_index) # affect!
    node_pars, edge_pars = integrator.p
    node_pars[event_index] .= [first(node_pars[event_index]), true]
    println("enter FRT at ", integrator.t, " @ ", event_index)
end
function ExitLocalErrorState(integrator, event_index) # affect!
    node_pars, edge_pars = integrator.p
    node_pars[event_index] .= [first(node_pars[event_index]), false]
    println("exit FRT at ", integrator.t, " @ ", event_index)
end

function LocalVoltageLimitDetection(UsollMV, vlim, volt_idxs)
    vl = LocalVoltageLimitCondition(UsollMV, vlim, volt_idxs)
    cb_affect = EnterLocalErrorState # ok
    cb_neg_affect = ExitLocalErrorState # ok
    return VectorContinuousCallback(vl, cb_affect, cb_neg_affect, length(volt_idxs)-1)
end


############### callback for voltage_limit_achieved #############

struct VoltageLimitCondition
    UsollMV
    vlim
    volt_idxs
    limit
end
function (vla::VoltageLimitCondition)(u, t, integrator) # condition
    # vla.UsollMV = 20. # kV
    # vla.vmax = 1.05 * UsollMV # 5% margin
    # vla.vmin = 0.95 * UsollMV
    # voltage_limit_achieved = 0

    # get voltage amplitudes (without slack bus on HV side)
    v_abs = [norm(u[v_idx]) for v_idx in vla.volt_idxs[2:end]]
    # println("max(v_abs) ", maximum(v_abs))
    #v_abs = [abs.(u[first.(v_idx)] .+ 1im .* u[last.(v_idx)]) for v_idx in vla.volt_idxs[2:end]]
    if vla.limit == :lower
        return minimum(v_abs) - vla.vlim
    elseif vla.limit == :upper
        return vla.vlim - maximum(v_abs)
    end
    #) || (minimum(v_abs) < vla.vmin) # we can always set to 1, -1 is not needed
    # upcrossing and downcrossing
end

function InitializeCallback!(c, u, t, integrator)
    #state = c.condition(u, t, integrator)
    #integrator.p = state == 0 ? 0 : 1
    # DO NOT CHANGE
    integrator.p = false
end

function StopErrorIntegrator!(integrator) # affect!
    node_pars, edge_pars = integrator.p
    for np in node_pars
        np .= [true, last(np)]
    end
    #integrator.p = true
    println("stop error integrator at ", integrator.t)
      # voltage_limit_achieved =
      # integrator.f.f.controller(integrator.u, integrator.t, integrator.f.f.volt_idxs, integrator.f.f.slack_connect, integrator.f.f.P_ref,  integrator.f.f.Q_ref)[3] = 1 # where is this in integrator? I can not find it in sol.prob.p
end

function RestartErrorIntegrator!(integrator) # affect!
    node_pars, edge_pars = integrator.p
    for np in node_pars
        np .= [false, last(np)]
    end
    #integrator.p = false
    println("restart error integrator at ", integrator.t)
      # voltage_limit_achieved =
      # integrator.f.f.controller(integrator.u, integrator.t, integrator.f.f.volt_idxs, integrator.f.f.slack_connect, integrator.f.f.P_ref,  integrator.f.f.Q_ref)[3] = 0 # where is this in integrator? I can not find it in sol.prob.p
end

# @Paul TODO Can you check this again? The integrator is stopped if the voltage achieves the limits and started within
# It looks the other way round here?!
# I adapted VoltageLimitDetection2(UsollMV, vlim_in, vlim_out, volt_idxs) which I think is correct
# I can not test it because neither of them is being used!


function VoltageLimitDetection(UsollMV, vlim_in, vlim_out, volt_idxs)
    # into error at upper bound
    vl_upper_in = VoltageLimitCondition(UsollMV, vlim_in*UsollMV, volt_idxs, :upper)
    # into error at lower bound
    vl_lower_in = VoltageLimitCondition(UsollMV, vlim_in*UsollMV, volt_idxs, :lower)
    # out of error at upper bound
    vl_upper_out = VoltageLimitCondition(UsollMV, vlim_out*UsollMV, volt_idxs, :upper)
    # out of error at lower bound
    vl_lower_out = VoltageLimitCondition(UsollMV, vlim_out*UsollMV, volt_idxs, :lower)
    return CallbackSet(
        ContinuousCallback(
            vl_lower_in,
            RestartErrorIntegrator!,
            nothing,  # needed, otherwise it is activated at up AND downcrossing
            #initialize = InitializeCallback!,
            rootfind=true,
            save_positions=(true,true),
              interp_points=10,
              abstol=1e-12,
              reltol=0,
              idxs=nothing
        ),
        ContinuousCallback(
            vl_upper_in,
            RestartErrorIntegrator!,
            nothing,
            #initialize = InitializeCallback!,
            rootfind=true,
              save_positions=(true,true),
              interp_points=10,
              abstol=1e-12,
              reltol=0,
              idxs=nothing
        ),
        ContinuousCallback(
            vl_lower_out,
            nothing,
            StopErrorIntegrator!,
            #initialize = InitializeCallback!,
            rootfind=true,
            save_positions=(true,true),
              interp_points=10,
              abstol=1e-12,
              reltol=0,
              idxs=nothing
        ),
        ContinuousCallback(
            vl_upper_out,
            nothing,
            StopErrorIntegrator!, # needed, otherwise it is activated at up AND downcrossing
            #initialize = InitializeCallback!,
            rootfind=true,
             save_positions=(true,true),
              interp_points=10,
              abstol=1e-12,
              reltol=0,
              idxs=nothing
        ),
    )
end
