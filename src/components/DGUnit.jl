#using LinearAlgebra
#using NetworkDynamics
import Base: @__doc__
using PowerDynamics: @DynamicNode
import PowerDynamics: construct_vertex, dimension, symbolsof, showdefinition, AbstractNode

@doc doc"""
```Julia
DGUnit(;I_r)
```
Pr_l_max = 2.0 # Holms Model: 10e12 / base_power
Rr_l_min = -2.0 # Holms Model: -10e12 / base_power
T = 10.0 # unit: s
K = 1.0 # unit: [y]/[u] = [P]/[P] = 1
Pmax = 2. # in this Limiter it is set to 200 MW in Holms Model!!
Pmin = -2. # in this limiter it is set to 0MW in Holms Model!
Qmax = 2. # in this Limiter it is set to 60 MW
Qmin = -2. # in this limiter it is set to -60 MW
Vac_ref = 1 * 20 * sqrt(2 / 3) # phase-ground = phase-phase * sqrt(2/3) in kV
V_dead = 0.1 * Vac_ref
k_FRT = 2.0 # p.u. or A/V or kA/kV stays the same as long as imax and Vac, Vac_ref and V_dead are all p.u or A,V or kA/kV
imax = 50 / sqrt(3) / base_current # assumption: 1 MVA nominal power of DGunit at 20kV and SB = 25 MVA
K_pll = 0.1 * base_voltage # theta in rad, Kpll before in rad/(Vs) now * 1000 V/kV --> rad/(kVs)
Y_shunt = 0.4 # shunt admittance for e.g. node short circuit
"""

begin
    @DynamicNode DGUnit(
        P_bkgrnd,
        Q_bkgrnd,
        Pr_l_min,
        Pr_l_max,
        Qr_l_min,
        Qr_l_max,
        Pr_del_min,
        Pr_del_max,
        Qr_del_min,
        Qr_del_max,
        T_PT1,
        K_PT1,
        Vac_ref,
        V_dead,
        k_FRT,
        imax,
        K_pll,
        Y_shunt,
    ) begin
        Diagonal([0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        #MassMatrix(m_u = false,m_int = [true, true, true, true, true])
    end begin end [
        [theta, dtheta],
        [Pr_del, dPr_del],
        [Qr_del, dQr_del],
        [Pr_l, dPr_l],
        [Qr_l, dQr_l]
    ] begin

        #= Obtaining the controller inputs =#
        ΔP, ΔQ, global_limit, local_limit = p

        # freeze integration when global_limit == 1
        Psat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit
        Qsat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit

        ############ switch upon global_limit ############

        if global_limit
            dPr_l = 0.
            dQr_l = 0.
        else
            dPr_l = SmoothTriggeredLimiter(
                -ΔP,
                Pr_l,
                Pr_l_min,
                Pr_l_max,
                Psat) #ΔP

             dQr_l = SmoothTriggeredLimiter(
                -ΔQ,
                Qr_l,
                Qr_l_min,
                Qr_l_max,
                Qsat) #ΔQ
        end

        P_LC_ref = Pr_l * 0.5  # integrator slow down only
        Q_LC_ref = Qr_l * 0.5 # integrator slow down only

        dPr_del = PT1(Pr_del, P_LC_ref, K_PT1, T_PT1) # PT1: T * dy + y = K * u --> dy = 1/T * (K * u - y) --> dy = PT1(y, u, K, T)
        dQr_del = PT1(Qr_del, Q_LC_ref, K_PT1, T_PT1)

        Pr_g = SmoothClamp(Pr_del, Pr_del_max, Pr_del_min)  # SmoothClamp(Pr_del, ymax, ymin)
        Qr_g = SmoothClamp(Qr_del, Qr_del_max, Qr_del_min)  # SmoothClamp(Qr_del, ymax, ymin)

        ########### switch upon local limit ############
        Va_amp = abs(u)
        id =  Pr_g / Va_amp # 2.0 / 3.0 # kA
        iq = - Qr_g / Va_amp # 2.0 / 3.0 * # kA

        DeltaV = Vac_ref - Va_amp
        #local_limit = abs(DeltaV) > V_dead ? true : false
        if local_limit #fault
            println("DeltaV", DeltaV)
            iq_plus = 0.
            if DeltaV > V_dead
                iq_plus = - k_FRT * (DeltaV - V_dead)
            elseif DeltaV < -V_dead
                iq_plus = - k_FRT * (DeltaV + V_dead)
            end
            iq_r = SmoothClamp(iq_plus, imax)
            id_r = SmoothClamp(
                iq,
                sqrt(abs(imax^2 - iq_r^2))
            )
        else #normal
            id_r = SmoothClamp(id, imax)
            iq_r = SmoothClamp(
                iq,
                sqrt(abs(imax^2 - id_r^2))
            )
        end

        # PLL
        dtheta = PLL(K_pll, u, theta)
        I_r = exp(1im * theta) * complex(id_r, iq_r)
        # println("I_r ", I_r)
        I_bg = conj(complex(P_bkgrnd, Q_bkgrnd)) / conj(u)
        # println("I_bg ", I_bg)
        # println("u ", u)
        # println("P_bkgrnd ", P_bkgrnd)
        # println("Q_bkgrnd ", Q_bkgrnd)
        I_shunt = - Y_shunt(t) * (0. - u) # total current and static line give these signs when thinking of the shunt as a line with source at node and destination at ground
        # println("I_shunt", I_shunt)
        # println("I_sum", I_r + I_bg + I_shunt)

        #println(I_bg)
        ### CSI ####
        du = (i + I_shunt) - (I_r + I_bg) # mass matrix = 0
        end
end




###############################
##### helping functions ######
###############################

function PT1(var, inp, K, T)  # dy = 1/T * (K * u - y)
    if T == 0
        error("T = 0 is not a valid input")
    else
        1.0 / T * (K * inp - var) # delay of 0.01 missing
    end
end

function smooth_cutoff(value, threshold; scale=10.)
    log(exp(scale * value) + exp.(scale * threshold)) / scale
end

function smooth_step(x, loc; order=10)
    (1 .+tanh.(order*(x-loc))) / 2.
end

# x = -2:.01:2
# plot(x, SmoothStep.(x, 1, -1))
# plot!(x, SmoothClamp.(x, 1, -1))
function SmoothStep(x, upper, lower; order=20)
    norm_fac =  1. / maximum(abs, (lower, upper))
    (1. - smooth_step.(x * norm_fac, upper * norm_fac; order=order)) .* (1. .- smooth_step.(-x * norm_fac, -lower * norm_fac; order=order))
end
SmoothStep(x, threshold) = SmoothStep(x, threshold, -threshold; order=10)

function Step(x, upper, lower)
    return lower < x < upper ? 1 : 0
end
 Step(x, threshold) =  Step(x, threshold, -threshold)

function SmoothClamp(value, upper, lower; scale=20.)
    if upper == lower
        return upper
    else
        norm_fac =  1. / maximum(abs, (lower, upper))
        if value > 0
            return -smooth_cutoff(-value * norm_fac, -upper * norm_fac; scale=scale) / norm_fac
        else
            return smooth_cutoff(value * norm_fac, lower * norm_fac; scale=scale) / norm_fac
        end
    end
end

SmoothClamp(value, threshold) = SmoothClamp(value, threshold, -threshold; scale = 20.)
SmoothClamp(value, threshold; scale=100.) = SmoothClamp(value, threshold, -threshold; scale = scale)


function SmoothFRTDetectionPlusLimit(err, Vac, Vac_ref, V_dead, k_FRT, id, iq, imax)
    DeltaV = Vac_ref - Vac
    #err = abs(DeltaV) > V_dead ? true : false
    iq_plus = 0.
        # THIS SHOULD STAY
    if err
        if DeltaV > V_dead
            iq_plus = - k_FRT * (DeltaV - V_dead)
        elseif DeltaV < -V_dead
            iq_plus = - k_FRT * (DeltaV + V_dead)
        end
        iq_r_fault = SmoothClamp(iq_plus, imax)
        id_r_fault = SmoothClamp(
            iq,
            sqrt(abs(imax^2 - iq_r_fault^2))
        )
        return  (id_r_fault, iq_r_fault)
    else
        iq_plus = 0.
        id_r_normal = SmoothClamp(id, imax)
        iq_r_normal = SmoothClamp(
            iq,
            sqrt(abs(imax^2 - id_r_normal^2))
        )
        return (id_r_normal, iq_r_normal)
    end
end

function SmoothFRTDetectionPlusLimit_added_iq(Vac, Vac_ref, V_dead, k_FRT, id, iq, imax)
    DeltaV = Vac_ref - Vac
    err = abs(DeltaV) > V_dead ? true : false

        # THIS SHOULD STAY
    if err
        if DeltaV > V_dead
            iq_plus = - k_FRT * (DeltaV - V_dead)
        elseif DeltaV < -V_dead
            iq_plus = - k_FRT * (DeltaV + V_dead)
        end
        iq_r_fault = SmoothClamp(iq + iq_plus, imax)
        id_r_fault = SmoothClamp(iq +
            iq_plus,
            sqrt(abs(imax^2 - iq_r_fault^2))
        )
        return  (id_r_fault, iq_r_fault)
    else
        iq_plus = 0
        id_r_normal = SmoothClamp(id, imax)
        iq_r_normal = SmoothClamp(
            iq,
            sqrt(abs(imax^2 - id_r_normal^2))
        )
        return (id_r_normal, iq_r_normal)
    end
end


function SmoothTriggeredLimiter(var, trigger, trigger_min, trigger_max, var_sat)
    if trigger_max < trigger_min
        error("upper limit $trigger_max, of var $var needs to be larger than lower limit $trigger_min")
    end
    if var_sat != 0
        error("TODO implement nonzero saturation value. Input: $var_sat")
    end
        return var * SmoothStep(trigger, trigger_max, trigger_min)
end




##############################################################
######### non smooth Limiter and old FRT #########################
#############################################





function TriggeredLimiter(var, trigger, trigger_min, trigger_max, var_sat, limit_achieved)
    if trigger_max < trigger_min
        error("upper limit $trigger_max, of var $var needs to be larger than lower limit $trigger_min")
    end
    out = var
    if trigger > trigger_max
        out = var_sat # 0
    elseif trigger < trigger_min
        out = var_sat # 0
    elseif limit_achieved != 0.
        out = var_sat # 0
    else
        out = var
    end
    out
end



function Limiter(var::Number, max, min)
    if max < min
        error("upper limit $max, of var $var needs to be larger than lower limit $min")
    end
    if var >= max
        var = max
        #println("max limit of $var achieved")
    elseif var <= min
        var = min
        #println("min limit of $var achieved")
    end
    var
end

function PT1(var, inp, K, T)  # dy = 1/T * (K * u - y)
    if T == 0
        error("T = 0 is not a valid input")
    else
        1.0 / T * (K * inp - var) # delay of 0.01 missing
    end
end

function PLL(K_pll, u, theta)
    #vq = (imag(u) * cos(theta) - real(u) * sin(theta)) # in simulink imag(u) and real(u) are calculted from V_amp and V_phase
    # this is 0 when angle(u)==theta
    K_pll * imag(u * exp(-1im * theta))
end


""" FRTDetection(Vac, Vac_ref, V_dead, k_FRT)
calculates if the voltage deviation is too strong and FRT is needed.
"""

# Vac_ref = 1.0
# V_dead = 0.1
# k_FRT = 2
function FRTDetection(Vac, Vac_ref, V_dead, k_FRT)

    DeltaV = Vac_ref - Vac
    #print(DeltaV)

    #if Vac > (Vac_ref + V_dead) || Vac < Vac_ref - V_dead
    err = abs(DeltaV) > V_dead ? 1 : 0

    # THIS SHOULD STAY
    if DeltaV > V_dead
        iq_plus = - k_FRT * (DeltaV - V_dead)
    elseif DeltaV < -V_dead
        iq_plus = - k_FRT * (DeltaV + V_dead)
    else
        iq_plus = 0
    end
    (err, iq_plus) # iq_plus != 0 also if error 0? Yes, in FRTLimit this is rechanged
end

""" FRTLimit(id,error,iq_plus,iq,imax)
calculates iq_r and id_r.
"""

function FRTLimit(id, err, iq_plus, iq, imax)
    # return id, iq # debugging
    id_r_normal = Limiter(id, imax, -imax)
    iq_r_normal = Limiter(
        iq,
        sqrt(abs(imax^2 - id_r_normal^2)),
        -sqrt(abs(imax^2 - id_r_normal^2)),
    )
    iq_r_fault = Limiter(iq_plus, imax, -imax)
    id_r_fault = Limiter(
        iq,
        sqrt(abs(imax^2 - iq_r_fault^2)),
        -sqrt(abs(imax^2 - iq_r_fault^2)),
    )


    if err == 1
        id_r = id_r_fault
        iq_r = iq_r_fault
    else
        id_r = id_r_normal
        iq_r = iq_r_normal
    end

    (id_r, iq_r)
end


function SmoothFRTLimit(id, err, iq_plus, iq, imax)

    # return id, iq # debugging
    # all clamps are symmetric
    id_r_normal = SmoothClamp(id, imax)
    iq_r_normal = SmoothClamp(
        iq,
        sqrt(abs(imax^2 - id_r_normal^2))
    )
    iq_r_fault = SmoothClamp(iq_plus, imax)
    id_r_fault = SmoothClamp(
        iq,
        sqrt(abs(imax^2 - iq_r_fault^2))
    )
    if err == 1
        return  (id_r_fault, iq_r_fault)
    else
        return (id_r_normal, iq_r_normal)
    end
end




##############################################################
##############################################################
### WRAPPER functions
##############################################################
###############################################################


""" RefSignal(SimulationTime, StepSize, ReferenceChanges)
calculates a Reference signal over the full simulation time (tested by Lia)
"""


function RefSignal(SimulationTime, StepSize, ReferenceChanges)
        #function [ ReferenceSignal ] = RefSignal( SimulationTime, StepSize, ReferenceChanges)

    ## Data input
    # Set simulation parameters
    SimulationTime = SimulationTime   # simulation time in seconds
    StepSize = StepSize      # step size in seconds


    # Define the values and timesteps of reference value changes
    # ReferenceChanges = [
    #     2   25*1e6
    #     4.2 10*1e6
    #     6   31*1e6
    #     9   15*1e6
    #     ];


    ## Initialization
    ReferenceChanges = [ReferenceChanges; [SimulationTime ReferenceChanges[end, 2]]]

    Time = collect(0:StepSize:SimulationTime)       # vector of time steps
    SizeTimeSteps = size(Time)                # size of vector of time steps
    ReferenceSignal = zeros(length(Time), 2) # added by Lia
    ReferenceSignal[:, 1] = Time
    ReferenceSignal[:, 2] = ones(SizeTimeSteps)


    ## Create Signal
    for it = 1:size(ReferenceChanges, 1)-1 # Loop for each change of reference signal

        IntervalStart = ReferenceChanges[it, 1] # ok, if this is a vector?
        IntervalEnd = ReferenceChanges[it+1, 1]

        StartIndex = findall(ReferenceSignal[:, 1] .<= IntervalStart)
        StartIndex = StartIndex[end]
        EndIndex = findall(ReferenceSignal[:, 1] .<= IntervalEnd)
        EndIndex = EndIndex[end]

        IntervalIndex = StartIndex:EndIndex

        ReferenceSignal[IntervalIndex, 2] .= ReferenceChanges[it, 2]
        ReferenceSignal
    end
end



""" VoltageLimitAchieve(V_amp, limit_achieve_old)
V_amp is here still the vector of amplitudes of the voltages at all DG units. The function checks if the voltage limit has been achieved at any DG unit. Will be used here to freeze the integrator!
    Unclear if the limit_achieve_old is needed as input.
"""


function VoltageLimitAchieve(V_amp, limit_achieve_old)
    # is limit_achieve_old necessary as input? if condition does the same (1.04/1.05 or 0.96/0.95??)
    #    function voltage_limit_achive   = fcn(V_amp, limit_achive_old)
    ##codegen

    u_max = maximum(V_amp)

    u_min = minimum(V_amp)
    # Voltage smaller than minimum
    if limit_achive_old == -1

        # Voltage smaller than minimum
        if u_min < 0.96

            voltage_limit_achive = -1


            # Voltage bigger than maximum
        elseif u_max > 1.05

            voltage_limit_achive = 1

        else

            # Voltage within the voltage boundaries
            voltage_limit_achive = 0


        end

    # Voltage bigger than maximum
    elseif limit_achive_old == 1

        # Voltage smaller than minimum
        if u_min < 0.95

            voltage_limit_achive = -1


            # Voltage bigger than maximum
        elseif u_max > 1.04

            voltage_limit_achive = 1

        else

            # Voltage within the voltage boundaries
            voltage_limit_achive = 0


        end


    else
        # Voltage smaller than minimum
        if u_min < 0.95

            voltage_limit_achive = -1


            # Voltage bigger than maximum
        elseif u_max > 1.05

            voltage_limit_achive = 1

        else

            # Voltage within the voltage boundaries
            voltage_limit_achive = 0


        end
    end
    voltage_limit_achive
end



function LimitAchieve(max, u, min)
        # function limit = fcn(max,u,min)
    ##codegen
    #P_DG_ref<P_DG_min
    if u < min
        limit = -1
    #P_DG_ref>P_DG_max
    elseif u > max
        limit = 1
    #P_DG_ref is ok
    else
        limit = 0
    end
    limit
end




# function ErrorReset(e,power_lim_reach,quantyn,dPQ)
#
#     e_out = 0;
#     reset = 0;
#
#     # P_ref < P_DG_min => P_lim -1 (positive e_P_out still allowed)
#     # P_ref > P_DG_max => P_lim 1  (negative e_P_out still allowed)
#     if power_lim_reach == -1
#         if e <= 0; e_out = 0; end
#     elseif power_lim_reach == 1
#         if e >= 0; e_out = 0; end
#     else
#         if quantyn == 0
#             e_out = e;
#         else
#             if abs(e) < dPQ
#                 e_out = 0;
#             else
#                 e_out = e;
#             end
#         end
#     end
#     (e_out,reset)
# end


##########################################################################################
########################## DGUnit_simple #################################################
##########################################################################################

begin
    @DynamicNode DGUnit_simple(
        P_bkgrnd,
        Q_bkgrnd,
        Pr_l_min,
        Pr_l_max,
        Qr_l_min,
        Qr_l_max,
        Pr_del_min,
        Pr_del_max,
        Qr_del_min,
        Qr_del_max,
        T_PT1,
        K_PT1,
        Vac_ref,
        V_dead,
        k_FRT,
        imax,
        K_pll,
        Y_shunt,
    ) begin
        Diagonal([0.0, 0.0, 1.0, 1.0])
        #MassMatrix(m_u = false,m_int = [true, true, true, true, true])
    end begin end [
    #    [theta, dtheta],
    #    [Pr_del, dPr_del],
    #    [Qr_del, dQr_del],
        [Pr_l, dPr_l],
        [Qr_l, dQr_l]
    ] begin

        #= Obtaining the controller inputs =#
        ΔP, ΔQ, voltage_limit_achieved, local_limit = p

        # freeze integration when voltage_limit_achieved == 1
        Psat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit
        Qsat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit

        # return when voltage_limit_achieved
        if voltage_limit_achieved
            dPr_l = 0.
            dQr_l = 0.
        else
            dPr_l = SmoothTriggeredLimiter(
                ΔP,
                Pr_l,
                Pr_l_min,
                Pr_l_max,
                Psat#,
                #voltage_limit_achieved,
            ) #ΔP

             dQr_l = SmoothTriggeredLimiter(
                ΔQ,
                Qr_l,
                Pr_l_min,
                Pr_l_max,
                Qsat#,
                #voltage_limit_achieved,
            ) #ΔQ
        end

        P_LC_ref = Pr_l * -0.5  # integrator slow down only
        Q_LC_ref = Qr_l * -0.5 # integrator slow down only

        #dPr_del = PT1(Pr_del, P_LC_ref, K_PT1, T_PT1) # PT1: T * dy + y = K * u --> dy = 1/T * (K * u - y) --> dy = PT1(y, u, K, T)
        #dQr_del = PT1(Qr_del, Q_LC_ref, K_PT1, T_PT1)

        #Pr_g = SmoothClamp(Pr_del, Pr_del_max, Pr_del_min)  # SmoothClamp(Pr_del, ymax, ymin)
        #Qr_g = SmoothClamp(Qr_del, Qr_del_max, Qr_del_min)  # SmoothClamp(Qr_del, ymax, ymin)

        Pr_g = P_LC_ref
        Qr_g = Q_LC_ref

        Va_amp = abs(u)

        id =  Pr_g / Va_amp # 2.0 / 3.0 # kA
        iq =  - Qr_g / Va_amp # 2.0 / 3.0 * # kA

        # (id_r, iq_r) = SmoothFRTDetectionPlusLimit(#_added_iq(
        #     Va_amp,
        #     Vac_ref,
        #     V_dead,
        #     k_FRT,
        #     id,
        #     iq,
        #     imax,
        # )
        # no FRT
        id_r = id
        iq_r = iq

        # PLL
    #    dtheta = PLL(K_pll, u, theta)
        I_r = complex(id_r, iq_r) # exp(1im * theta)
        I_bg = conj(complex(P_bkgrnd, Q_bkgrnd)) / conj(u)
        I_shunt = Y_shunt(t) * u
        #println(I_bg)
        ### CSI ####
        du = i - (I_r + I_bg + I_shunt) # mass matrix = 0
        end
end




##########################################################################################
########################## DGUnit_simple_FRT #################################################
##########################################################################################

begin
    @DynamicNode DGUnit_simple_FRT(
        P_bkgrnd,
        Q_bkgrnd,
        Pr_l_min,
        Pr_l_max,
        Qr_l_min,
        Qr_l_max,
        Pr_del_min,
        Pr_del_max,
        Qr_del_min,
        Qr_del_max,
        T_PT1,
        K_PT1,
        Vac_ref,
        V_dead,
        k_FRT,
        imax,
        K_pll,
        Y_shunt,
    ) begin
        Diagonal([0.0, 0.0, 1.0, 1.0])
        #MassMatrix(m_u = false,m_int = [true, true, true, true, true])
    end begin end [
    #    [theta, dtheta],
    #    [Pr_del, dPr_del],
    #    [Qr_del, dQr_del],
        [Pr_l, dPr_l],
        [Qr_l, dQr_l]
    ] begin

        #= Obtaining the controller inputs =#
        ΔP, ΔQ, voltage_limit_achieved, local_limit = p

        # freeze integration when voltage_limit_achieved == 1
        Psat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit
        Qsat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit

        # return when voltage_limit_achieved
        if voltage_limit_achieved
            dPr_l = 0.
            dQr_l = 0.
        else
            dPr_l = SmoothTriggeredLimiter(
                ΔP,
                Pr_l,
                Pr_l_min,
                Pr_l_max,
                Psat#,
                #voltage_limit_achieved,
            ) #ΔP

             dQr_l = SmoothTriggeredLimiter(
                ΔQ,
                Qr_l,
                Pr_l_min,
                Pr_l_max,
                Qsat#,
                #voltage_limit_achieved,
            ) #ΔQ
        end

        P_LC_ref = Pr_l * -0.5  # integrator slow down only
        Q_LC_ref = Qr_l * -0.5 # integrator slow down only

        #dPr_del = PT1(Pr_del, P_LC_ref, K_PT1, T_PT1) # PT1: T * dy + y = K * u --> dy = 1/T * (K * u - y) --> dy = PT1(y, u, K, T)
        #dQr_del = PT1(Qr_del, Q_LC_ref, K_PT1, T_PT1)

        #Pr_g = SmoothClamp(Pr_del, Pr_del_max, Pr_del_min)  # SmoothClamp(Pr_del, ymax, ymin)
        #Qr_g = SmoothClamp(Qr_del, Qr_del_max, Qr_del_min)  # SmoothClamp(Qr_del, ymax, ymin)

        Pr_g = P_LC_ref
        Qr_g = Q_LC_ref

        Va_amp = abs(u)

        id =  Pr_g / Va_amp # 2.0 / 3.0 # kA
        iq =  - Qr_g / Va_amp # 2.0 / 3.0 * # kA

        DeltaV = Vac_ref - Va_amp
        if local_limit #fault
            iq_plus = 0.
            if DeltaV > V_dead
                iq_plus = - k_FRT * (DeltaV - V_dead)
            elseif DeltaV < -V_dead
                iq_plus = - k_FRT * (DeltaV + V_dead)
            end
            iq_r = SmoothClamp(iq_plus, imax)
            id_r = SmoothClamp(
                iq,
                sqrt(abs(imax^2 - iq_r^2))
            )
        else #normal
            id_r = SmoothClamp(id, imax)
            iq_r = SmoothClamp(
                iq,
                sqrt(abs(imax^2 - id_r^2))
            )
        end

        # PLL
    #    dtheta = PLL(K_pll, u, theta)
        I_r = complex(id_r, iq_r) # exp(1im * theta)
        I_bg = conj(complex(P_bkgrnd, Q_bkgrnd)) / conj(u)
        I_shunt = Y_shunt(t) * u
        #println(I_bg)
        ### CSI ####
        du = i - (I_r + I_bg + I_shunt) # mass matrix = 0
        end
end



##########################################################################################
########################## DGUnit_simple_PLL #################################################
##########################################################################################

begin
    @DynamicNode DGUnit_simple_PLL(
        P_bkgrnd,
        Q_bkgrnd,
        Pr_l_min,
        Pr_l_max,
        Qr_l_min,
        Qr_l_max,
        Pr_del_min,
        Pr_del_max,
        Qr_del_min,
        Qr_del_max,
        T_PT1,
        K_PT1,
        Vac_ref,
        V_dead,
        k_FRT,
        imax,
        K_pll,
        Y_shunt,
    ) begin
        Diagonal([0.0, 0.0, 1.0, 1.0, 1.0])
        #MassMatrix(m_u = false,m_int = [true, true, true, true, true])
    end begin end [
        [theta, dtheta],
    #    [Pr_del, dPr_del],
    #    [Qr_del, dQr_del],
        [Pr_l, dPr_l],
        [Qr_l, dQr_l]
    ] begin

        #= Obtaining the controller inputs =#
        ΔP, ΔQ, voltage_limit_achieved, local_limit = p

        # freeze integration when voltage_limit_achieved == 1
        Psat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit
        Qsat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit

        # return when voltage_limit_achieved
        if voltage_limit_achieved
            dPr_l = 0.
            dQr_l = 0.
        else
            dPr_l = SmoothTriggeredLimiter(
                ΔP,
                Pr_l,
                Pr_l_min,
                Pr_l_max,
                Psat#,
                #voltage_limit_achieved,
            ) #ΔP

             dQr_l = SmoothTriggeredLimiter(
                ΔQ,
                Qr_l,
                Pr_l_min,
                Pr_l_max,
                Qsat#,
                #voltage_limit_achieved,
            ) #ΔQ
        end

        P_LC_ref = Pr_l * -0.5  # integrator slow down only
        Q_LC_ref = Qr_l * -0.5 # integrator slow down only

        #dPr_del = PT1(Pr_del, P_LC_ref, K_PT1, T_PT1) # PT1: T * dy + y = K * u --> dy = 1/T * (K * u - y) --> dy = PT1(y, u, K, T)
        #dQr_del = PT1(Qr_del, Q_LC_ref, K_PT1, T_PT1)

        #Pr_g = SmoothClamp(Pr_del, Pr_del_max, Pr_del_min)  # SmoothClamp(Pr_del, ymax, ymin)
        #Qr_g = SmoothClamp(Qr_del, Qr_del_max, Qr_del_min)  # SmoothClamp(Qr_del, ymax, ymin)

        Pr_g = P_LC_ref
        Qr_g = Q_LC_ref

        Va_amp = abs(u)

        id =  Pr_g / Va_amp # 2.0 / 3.0 # kA
        iq =  - Qr_g / Va_amp # 2.0 / 3.0 * # kA

        id_r = id
        iq_r = iq

        # PLL
       dtheta = PLL(K_pll, u, theta)
        I_r = complex(id_r, iq_r) * exp(1im * theta)
        I_bg = conj(complex(P_bkgrnd, Q_bkgrnd)) / conj(u)
        I_shunt = Y_shunt(t) * u
        #println(I_bg)
        ### CSI ####
        du = i - (I_r + I_bg + I_shunt) # mass matrix = 0
        end
end


##########################################################################################
########################## DGUnit_nolimit_new #################################################
##########################################################################################

begin
    @DynamicNode DGUnit_nolimit_new(
        P_bkgrnd,
        Q_bkgrnd,
        Pr_l_min,
        Pr_l_max,
        Qr_l_min,
        Qr_l_max,
        Pr_del_min,
        Pr_del_max,
        Qr_del_min,
        Qr_del_max,
        T_PT1,
        K_PT1,
        Vac_ref,
        V_dead,
        k_FRT,
        imax,
        K_pll,
        Y_shunt,
    ) begin
        Diagonal([0.0, 0.0, 1.0, 1.0])
        #MassMatrix(m_u = false,m_int = [true, true, true, true, true])
    end begin end [
    #    [theta, dtheta],
    #    [Pr_del, dPr_del],
    #    [Qr_del, dQr_del],
        [Pr_l, dPr_l],
        [Qr_l, dQr_l]
    ] begin

        #= Obtaining the controller inputs =#
        ΔP, ΔQ, voltage_limit_achieved, local_limit = p

        # return when voltage_limit_achieved
        if voltage_limit_achieved
            dPr_l = 0.
            dQr_l = 0.
        else
            dPr_l = ΔP
            dQr_l = ΔQ
        end

        I_r = complex(Pr_l * -0.5 / abs(u),  - Qr_l * -0.5  / abs(u))
        I_bg = conj(complex(P_bkgrnd, Q_bkgrnd)) / conj(u)
        I_shunt = Y_shunt(t) * u
        #println(I_bg)
        ### CSI ####
        du = i - (I_r + I_bg + I_shunt) # mass matrix = 0
        end
end



##########################################################################################
########################## DGUnit_nolimit_FRT #################################################
##########################################################################################

begin
    @DynamicNode DGUnit_nolimit_FRT(
        P_bkgrnd,
        Q_bkgrnd,
        Pr_l_min,
        Pr_l_max,
        Qr_l_min,
        Qr_l_max,
        Pr_del_min,
        Pr_del_max,
        Qr_del_min,
        Qr_del_max,
        T_PT1,
        K_PT1,
        Vac_ref,
        V_dead,
        k_FRT,
        imax,
        K_pll,
        Y_shunt,
    ) begin
        Diagonal([0.0, 0.0, 1.0, 1.0])
        #MassMatrix(m_u = false,m_int = [true, true, true, true, true])
    end begin end [
    #    [theta, dtheta],
    #    [Pr_del, dPr_del],
    #    [Qr_del, dQr_del],
        [Pr_l, dPr_l],
        [Qr_l, dQr_l]
    ] begin

        #= Obtaining the controller inputs =#
        ΔP, ΔQ, voltage_limit_achieved, local_limit = p

        # return when voltage_limit_achieved
        if voltage_limit_achieved
            dPr_l = 0.
            dQr_l = 0.
        else
            dPr_l = ΔP
            dQr_l = ΔQ
        end

        id = Pr_l * -0.5 / abs(u)
        iq = - Qr_l * -0.5  / abs(u)

        DeltaV = Vac_ref - abs(u)
        if local_limit #fault
            iq_plus = 0.
            if DeltaV > V_dead
                iq_plus = - k_FRT * (DeltaV - V_dead)
            elseif DeltaV < -V_dead
                iq_plus = - k_FRT * (DeltaV + V_dead)
            end
            iq_r = SmoothClamp(iq_plus, imax)
            id_r = SmoothClamp(
                iq,
                sqrt(abs(imax^2 - iq_r^2))
            )
        else #normal
            id_r = SmoothClamp(id, imax)
            iq_r = SmoothClamp(
                iq,
                sqrt(abs(imax^2 - id_r^2))
            )
        end

        I_r = complex(id_r,  iq_r)
        I_bg = conj(complex(P_bkgrnd, Q_bkgrnd)) / conj(u)
        I_shunt = Y_shunt(t) * u
        #println(I_bg)
        ### CSI ####
        du = i - (I_r + I_bg + I_shunt) # mass matrix = 0
        end
end

##########################################################################################
########################## DGUnit_nolimit_PLL #################################################
##########################################################################################

begin
    @DynamicNode DGUnit_nolimit_PLL(
        P_bkgrnd,
        Q_bkgrnd,
        Pr_l_min,
        Pr_l_max,
        Qr_l_min,
        Qr_l_max,
        Pr_del_min,
        Pr_del_max,
        Qr_del_min,
        Qr_del_max,
        T_PT1,
        K_PT1,
        Vac_ref,
        V_dead,
        k_FRT,
        imax,
        K_pll,
        Y_shunt,
    ) begin
        Diagonal([0.0, 0.0, 1.0, 1.0, 1.0])
        #MassMatrix(m_u = false,m_int = [true, true, true, true, true])
    end begin end [
        [theta, dtheta],
    #    [Pr_del, dPr_del],
    #    [Qr_del, dQr_del],
        [Pr_l, dPr_l],
        [Qr_l, dQr_l]
    ] begin

        #= Obtaining the controller inputs =#
        ΔP, ΔQ, voltage_limit_achieved, local_limit = p

        # freeze integration when voltage_limit_achieved == 1
        Psat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit
        Qsat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit

        if voltage_limit_achieved
            dPr_l = 0.
            dQr_l = 0.
        else
            dPr_l = ΔP
            dQr_l = ΔQ
        end

        P_LC_ref = Pr_l * -0.5  # integrator slow down only
        Q_LC_ref = Qr_l * -0.5 # integrator slow down only

        #dPr_del = PT1(Pr_del, P_LC_ref, K_PT1, T_PT1) # PT1: T * dy + y = K * u --> dy = 1/T * (K * u - y) --> dy = PT1(y, u, K, T)
        #dQr_del = PT1(Qr_del, Q_LC_ref, K_PT1, T_PT1)

        #Pr_g = SmoothClamp(Pr_del, Pr_del_max, Pr_del_min)  # SmoothClamp(Pr_del, ymax, ymin)
        #Qr_g = SmoothClamp(Qr_del, Qr_del_max, Qr_del_min)  # SmoothClamp(Qr_del, ymax, ymin)

        Pr_g = P_LC_ref
        Qr_g = Q_LC_ref

        Va_amp = abs(u)

        id =  Pr_g / Va_amp # 2.0 / 3.0 # kA
        iq =  - Qr_g / Va_amp # 2.0 / 3.0 * # kA

        id_r = id
        iq_r = iq

        # PLL
       dtheta = PLL(K_pll, u, theta)
        I_r = complex(id_r, iq_r) * exp(1im * theta)
        I_bg = conj(complex(P_bkgrnd, Q_bkgrnd)) / conj(u)
        I_shunt = Y_shunt(t) * u
        #println(I_bg)
        ### CSI ####
        du = i - (I_r + I_bg + I_shunt) # mass matrix = 0
        end
end




##########################################################################################
########################## DGUnit_load #################################################
##########################################################################################

begin
    @DynamicNode DGUnit_load(
        P_bkgrnd,
        Q_bkgrnd,
        Pr_l_min,
        Pr_l_max,
        Qr_l_min,
        Qr_l_max,
        Pr_del_min,
        Pr_del_max,
        Qr_del_min,
        Qr_del_max,
        T_PT1,
        K_PT1,
        Vac_ref,
        V_dead,
        k_FRT,
        imax,
        K_pll,
        Y_shunt,
    ) begin
        Diagonal([0.0, 0.0])
        #MassMatrix(m_u = false,m_int = [true, true, true, true, true])
    end begin end [
    #    [theta, dtheta],
    #    [Pr_del, dPr_del],
    #    [Qr_del, dQr_del],
    #    [Pr_l, dPr_l],
    #    [Qr_l, dQr_l]
    ] begin

        #= Obtaining the controller inputs =#
        ΔP, ΔQ, voltage_limit_achieved, local_limit = p

        # return when voltage_limit_achieved
        # if voltage_limit_achieved
        #     dPr_l = 0.
        #     dQr_l = 0.
        # else
        #     dPr_l = ΔP
        #     dQr_l = ΔQ
        # end

        # id = Pr_l * -0.5 / abs(u)
        # iq = - Qr_l * -0.5  / abs(u)

        # (id_r, iq_r) = SmoothFRTDetectionPlusLimit(#_added_iq(
        #     abs(u),
        #     Vac_ref,
        #     V_dead,
        #     k_FRT,
        #     id,
        #     iq,
        #     imax,
        # )

        I_r = 0. #complex(id_r,  iq_r)
        I_bg = conj(complex(P_bkgrnd, Q_bkgrnd)) / conj(u)
        I_shunt = Y_shunt(t) * u
        #println(I_bg)
        ### CSI ####
        du = i - (I_r + I_bg + I_shunt) # mass matrix = 0
        end
end



##########################################################################################
################################ DG Unit nolimit OLD FROM FRANK ########################################
#########################################################################################


begin
    @DynamicNode DGUnit_nolimit(P_bkgrnd, Q_bkgrnd) begin
        Diagonal([0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0]) # MassMatrix(m_u = false,m_int = [false, false, false, true, true, true, true, true])
    end begin end [
        [theta, dtheta],
        [Pref, dPref],
        [Qref, dQref],
        [Pr_l, dPr_l],
        [Qr_l, dQr_l],
    ] begin

    #= Obtaining the controller inputs =#
        ΔP, ΔQ, voltage_limit_achieved, local_limit = p

    #= Begin of dynamics =#

    #= Integrating the error =#
        dPr_l = ΔP
        dQr_l = ΔQ

    #= PT1 =#
        T_inv = 0.1

        dPref = T_inv * (Pr_l - Pref)
        dQref = T_inv * (Qr_l - Qref)

    #= calculate current to inject in the corotating system =#
        Va_amp = abs(u) * sqrt(2/3)

        id_r = 2.0 / 3.0 * Pref / Va_amp
        iq_r = -2.0 / 3.0 * Qref / Va_amp

    #= PLL =#

        K_pll = 0.1
        vq = imag(u) * cos(theta) - real(u) * sin(theta) # in simulink imag(u) and real(u) are calculted from V_amp and V_phase
        dtheta = K_pll * vq

    #= CSI =#
        I_r =  exp(1.0im * theta) * complex(id_r, iq_r) + complex(P_bkgrnd, Q_bkgrnd) / u
        du = i - I_r # = 0. due to mass matrix
    end
end



##########################################################################################
################################ DG Unit dynamical ########################################
#########################################################################################
begin
    @DynamicNode DGUnit_dynamical(
        P_bkgrnd,
        Q_bkgrnd,
        Pr_l_min,
        Pr_l_max,
        Qr_l_min,
        Qr_l_max,
        Pr_del_min,
        Pr_del_max,
        Qr_del_min,
        Qr_del_max,
        T_PT1,
        K_PT1,
        Vac_ref,
        V_dead,
        k_FRT,
        imax,
        K_pll,
        Y_shunt,
    ) begin
        Diagonal([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        #MassMatrix(m_u = false,m_int = [true, true, true, true, true])
    end begin end [
        [theta, dtheta],
        [Pr_del, dPr_del],
        [Qr_del, dQr_del],
        [Pr_l, dPr_l],
        [Qr_l, dQr_l]
    ] begin
    #= Obtaining the controller inputs =#
    ΔP, ΔQ, global_limit, local_limit = p

    # freeze integration when global_limit == 1
    Psat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit
    Qsat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit

    ############ switch upon global_limit ############

    if global_limit
        dPr_l = 0.
        dQr_l = 0.
    else
        dPr_l = SmoothTriggeredLimiter(
            ΔP,
            Pr_l,
            Pr_l_min,
            Pr_l_max,
            Psat) #ΔP

         dQr_l = SmoothTriggeredLimiter(
            ΔQ,
            Qr_l,
            Qr_l_min,
            Qr_l_max,
            Qsat) #ΔQ
    end

    P_LC_ref = Pr_l * -0.5  # integrator slow down only
    Q_LC_ref = Qr_l * -0.5 # integrator slow down only

    dPr_del = PT1(Pr_del, P_LC_ref, K_PT1, T_PT1) # PT1: T * dy + y = K * u --> dy = 1/T * (K * u - y) --> dy = PT1(y, u, K, T)
    dQr_del = PT1(Qr_del, Q_LC_ref, K_PT1, T_PT1)

    Pr_g = SmoothClamp(Pr_del, Pr_del_max, Pr_del_min)  # SmoothClamp(Pr_del, ymax, ymin)
    Qr_g = SmoothClamp(Qr_del, Qr_del_max, Qr_del_min)  # SmoothClamp(Qr_del, ymax, ymin)

    ########### switch upon local limit ############
    Va_amp = abs(u)
    id =  Pr_g / Va_amp # 2.0 / 3.0 # kA
    iq =  - Qr_g / Va_amp # 2.0 / 3.0 * # kA

    DeltaV = Vac_ref - Va_amp
    if local_limit #fault
        iq_plus = 0.
        if DeltaV > V_dead
            iq_plus = - k_FRT * (DeltaV - V_dead)
        elseif DeltaV < -V_dead
            iq_plus = - k_FRT * (DeltaV + V_dead)
        end
        iq_r = SmoothClamp(iq_plus, imax)
        id_r = SmoothClamp(
            iq,
            sqrt(abs(imax^2 - iq_r^2))
        )
    else #normal
        id_r = SmoothClamp(id, imax)
        iq_r = SmoothClamp(
            iq,
            sqrt(abs(imax^2 - id_r^2))
        )
    end

    # PLL
    dtheta = PLL(K_pll, u, theta)
    I_r = exp(1im * theta) * complex(id_r, iq_r)
    I_bg = conj(complex(P_bkgrnd, Q_bkgrnd)) / conj(u)
    I_shunt = Y_shunt(t) * u

    ### CSI ####
    # introduce artifical capacity
    gain = 100.
    du = gain * (i -(I_r + I_bg + I_shunt))
    end
end


################################
########### DGUnit_without_PLL ###
###################################


begin
    @DynamicNode DGUnit_without_PLL(
        P_bkgrnd,
        Q_bkgrnd,
        Pr_l_min,
        Pr_l_max,
        Qr_l_min,
        Qr_l_max,
        Pr_del_min,
        Pr_del_max,
        Qr_del_min,
        Qr_del_max,
        T_PT1,
        K_PT1,
        Vac_ref,
        V_dead,
        k_FRT,
        imax,
        K_pll,
        Y_shunt,
    ) begin
        Diagonal([0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
        #MassMatrix(m_u = false,m_int = [true, true, true, true, true])
    end begin end [
    #    [theta, dtheta],
        [Pr_del, dPr_del],
        [Qr_del, dQr_del],
        [Pr_l, dPr_l],
        [Qr_l, dQr_l]
    ] begin

        #= Obtaining the controller inputs =#
        ΔP, ΔQ, global_limit, local_limit = p

        # freeze integration when global_limit == 1
        Psat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit
        Qsat = 0.0 # integrator set to zero if limits are triggered - should stay inside DGUnit

        ############ switch upon global_limit ############

        if global_limit
            dPr_l = 0.
            dQr_l = 0.
        else
            dPr_l = SmoothTriggeredLimiter(
                ΔP,
                Pr_l,
                Pr_l_min,
                Pr_l_max,
                Psat) #ΔP

             dQr_l = SmoothTriggeredLimiter(
                ΔQ,
                Qr_l,
                Qr_l_min,
                Qr_l_max,
                Qsat) #ΔQ
        end

        P_LC_ref = Pr_l * -0.5  # integrator slow down only
        Q_LC_ref = Qr_l * -0.5 # integrator slow down only

        dPr_del = PT1(Pr_del, P_LC_ref, K_PT1, T_PT1) # PT1: T * dy + y = K * u --> dy = 1/T * (K * u - y) --> dy = PT1(y, u, K, T)
        dQr_del = PT1(Qr_del, Q_LC_ref, K_PT1, T_PT1)

        Pr_g = Pr_del#SmoothClamp(Pr_del, Pr_del_max, Pr_del_min)  # SmoothClamp(Pr_del, ymax, ymin)
        Qr_g = Qr_del#SmoothClamp(Qr_del, Qr_del_max, Qr_del_min)  # SmoothClamp(Qr_del, ymax, ymin)

        ########### switch upon local limit ############
        Va_amp = abs(u)
        id =  Pr_g / Va_amp # 2.0 / 3.0 # kA
        iq =  - Qr_g / Va_amp # 2.0 / 3.0 * # kA

        DeltaV = Vac_ref - Va_amp
        #local_limit = abs(DeltaV) > V_dead ? true : false
        if local_limit #fault
            iq_plus = 0.
            if DeltaV > V_dead
                iq_plus = - k_FRT * (DeltaV - V_dead)
            elseif DeltaV < -V_dead
                iq_plus = - k_FRT * (DeltaV + V_dead)
            end
            iq_r = SmoothClamp(iq_plus, imax)
            id_r = SmoothClamp(
                iq,
                sqrt(abs(imax^2 - iq_r^2))
            )
        else #normal
            id_r = SmoothClamp(id, imax)
            iq_r = SmoothClamp(
                iq,
                sqrt(abs(imax^2 - id_r^2))
            )
        end

        # PLL
        # dtheta = PLL(K_pll, u, theta)
        I_r =  -complex(id_r, iq_r) #* exp(1im * theta)
        I_bg = conj(complex(P_bkgrnd, Q_bkgrnd)) / conj(u)
        I_shunt = Y_shunt(t) * u
        #println(I_bg)
        ### CSI ####
        du = i - (I_r - I_bg + I_shunt) # mass matrix = 0
        end
end
