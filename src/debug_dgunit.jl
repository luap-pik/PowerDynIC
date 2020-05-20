using Pkg
Pkg.activate(".")

using PowerDynamics
using OrdinaryDiffEq
using DifferentialEquations: DynamicSS

## set per unit MV
const base_power = 1E6 # 1MW
const base_voltage = 20E3 # 20kV
const base_current = base_power / base_voltage # 50A
const base_admittance = base_power / base_voltage^2 # 0.0025Ω^-1
const ω = 2 * π * 50.0 # 314.1593rad/s
# per unit HV
const base_voltage_HV = 110E3 # 110kV

## some additional functions that might develop into new PD.jl features
include("PDpatches.jl")
# control layer (needs ControlledPowerGrid from PDpatches)
include("components/control.jl")
# custom node type
include("components/DGUnit.jl")
include("components/OLTC.jl")
# TODO PD systemsize issue. base State type on GraphData??
include("components/RLLine.jl")
# load actual data
include("../example_cases/CIGRE/cigre_static.jl")


## define test cases

function testbed_chain(
    device,
    S,
    U = complex(1.0, 0.0),
    Y = complex(270.0, -380.0),
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
)
    busses = [SlackAlgebraic(; U = U), device, PQAlgebraic(; S = S)]
    lines = [
        StaticLine(; from = 1, to = 2, Y = Y),
        StaticLine(; from = 2, to = 3, Y = Y),
    ]
    pg = PowerGrid(busses, lines)
    pg, ADN(pg, typeof(device), P_ref, Q_ref)
end

function testbed_CIGRE(
    DG_unit_type;
    idx=1:3,
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
)
    busses_static, lines, T, elist, Zs, Yshs = CIGRE_static()
    busses = copy(busses_static)
    DG_locs = 2:12
    for i in DG_locs
        S_bkgrnd = zero(im)
        try
            S_bkgrnd = busses_static[i].S
        catch
            S_bkgrnd = complex(busses_static[i].P, busses_static[i].Q)
        end
        busses[i] = eval(Symbol(DG_unit_type))(;
            K_pll = 0.1base_voltage, #Hz/pu
            K_PT1 = 1.0, # unit: [y]/[u] = [P]/[P] = 1
            T_PT1 = 0.1, # unit: s
            S_pq = S_bkgrnd, # pu
            Y_n = 0.0,
        )
    end
    pg = PowerGrid(busses[idx], lines[1:length(idx)-1])
    pg, ADN(pg, DG_unit_type, P_ref, Q_ref)
end

function testbed(
    device1,
    device2,
    U = complex(1.0, 0.0),
    Y = complex(262, -375);
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
)
    T = OLTC(
        from = 1,
        to = 2,
        Uos = 110E3 / base_voltage_HV, # Bemessungsspannung Oberspannungsseite in kV
        Uus = 20E3 / base_voltage,# Bemessungsspannung Unterspannungsseite in kV
        k = 0, # Kennzahl der Schaltgruppe, Bsp: YD5-> 5
        ssp = 10.0, # Stufenschalterposition
        stufe_us = 0.625, # Stufung pro Stufenschalterstellung in %
        Sr = 25E6 / base_power, #Bemessungsscheinleistung in MVA
        uk = 12.0, # Kurzschlussspannung in %
        Pvk = 25E3 / base_power, # Kupferverluste in MW
        Pvl = 0.0, # Eisenverluste in kW
        iLeer = 0.0, # Leerlaufstrom in %
    )
    L = PiModelLine(
        from = 2,
        to = 3,
        y = Y,
        y_shunt_km = complex(0., 0.01),
        y_shunt_mk = complex(0., 0.01),
    )
    busses = [SlackAlgebraic(; U = U), device1, device2]
    lines = [T, L]
    pg = PowerGrid(busses, lines)
    pg, ADN(pg, typeof(device), P_ref, Q_ref)
end

function testbed_dynamicline(
    device,
    U = complex(1.0, 0.0),
    R = 0.5,
    L = 0.7;
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
)
    busses = [SlackAlgebraic(; U = U), device, device]
    lines = [
        RLLine(; from = 1, to = 2, R = R, L = L, ω0 = ω),
        RLLine(; from = 1, to = 2, R = R, L = L, ω0 = ω),
    ]
    pg = PowerGrid(busses, lines)
    pg, ADN(pg, typeof(device), P_ref, Q_ref)
end

function testbed_line(
    device,
    U = complex(1.0, 0.0),
    Y = complex(270.0, -380.0);
    P_ref = t -> -2,
    Q_ref = t -> 0.0,
)
    busses = [SlackAlgebraic(; U = U), device, device]
    lines = [
        StaticLine(; from = 1, to = 2, Y = Y),
        StaticLine(; from = 2, to = 3, Y = Y),
    ]
    pg = PowerGrid(busses, lines)
    pg, ADN(pg, typeof(device), P_ref, Q_ref)
end

# all values in p.u.
device = DGUnitPLLPQTracking(;
    K_pll = 0.1, #Hz/pu
    K_PT1 = 1.0, # unit: [y]/[u] = [P]/[P] = 1
    T_PT1 = 10.0, # unit: s
    S_pq = complex(-1.0, -0.25), # pu
    Y_n = 0.0,
)

## setup system

# define events
t_step = 10.0 # step in the power flow reference
t_fault = 1.0 # onset of node short circuit
t_duration = 0.1 # duration of node short circuit

P_ref(t) = t > t_step ? -1.0 : 0
Q_ref(t) = t > t_step ? -0.25 : 0

# edge data
# ode = rhs(pg)
# ode.f.graph_data.e

## find operating point

pg, cpg = testbed_CIGRE(DGUnitPLLPQTracking; idx=1:12, P_ref = t -> P_ref(0.0), Q_ref = t -> Q_ref(0.0))
ode = rhs(cpg)
state(vec) = State(pg, vec)
solution(vec) = PowerGridSolution(vec, pg)


op_guess = initial_guess(cpg)

#initial_error = zeros(2)
op_prob = ODEProblem(ode, op_guess, (0.0, t_step))#, initial_error)
_op = solve(
    SteadyStateProblem(op_prob),
    DynamicSS(Rodas5(); abstol = 1e-8, reltol = 1e-6, tspan=Inf)
)
op = _op |> state

## check dynamics

pg, cpg = testbed_CIGRE(DGUnitPLLPQTracking; idx=1:12, P_ref = P_ref, Q_ref = Q_ref)
ode = rhs(cpg)

# TODO: debug perturb
# TODO: why is there no change???
perturb = Perturbation(3, :θ, Dec(0.0))
u0 = perturb(op)

prob = ODEProblem(ode, u0.vec, (0.0, 30.0))
_sol = solve(
    prob,
    Rodas4(),
    d_discontinuities = [t_step, t_fault, t_fault + t_duration],
)
sol = _sol |> solution

## node short circuit

nsc = NodeShortCircuit(3, 10base_admittance, (t_fault, t_fault+t_duration)) # 10Ω

_sol = simulate_nsc(nsc, cpg, _op, (0., 30.))
sol = _sol |> solution

## plots

using Plots
plotly()

inc = (last(sol.dqsol.prob.tspan) - first(sol.dqsol.prob.tspan)) / 100.
t_array = range(first(sol.dqsol.prob.tspan); stop=last(sol.dqsol.prob.tspan), step=inc)


# TODO: bugs in solution object, e.g. (sol, 2, :v), (sol, 2, :φ)
# (sol, 2, :iabs) triggers call to rhs and hangs
# --> compare with indexing sol

plot(sol, :, :v, c=:gray, label=false)
plot!(sol, 2:3, :v)
hline!(op[2:3, :v], label = "op")
vline!([t_step,], label="step", c=:black)
vline!([t_fault, t_fault + t_duration], label="fault", c=:black)

plot(sol, :, :φ, c=:gray, label=false)
plot!(sol, 2:3, :φ)
plot!(sol, 2:3, :θ, c = "black", label = "PLL")
hline!(op[2:3, :φ], label = "op")
vline!([t_step,], label="step", c=:black)
vline!([t_fault, t_fault + t_duration], label="fault", c=:black)
#ylims!(-0.01, 0.01)

plot(sol, 2:3, :P_g)
hline!(op[2:3, :P_g], label = "op")
plot!(sol, 2:3, :Q_g)
hline!(op[2:3, :Q_g], label = "op")
vline!([t_step,], label="step", c=:black)
vline!([t_fault, t_fault + t_duration], label="fault", c=:black)


plot(sol, 2:3, :P_err)
plot!(sol, 2:3, :Q_err)

Δ = [
    cpg.controller(u, nothing, t) |> first
    for (t, u) in zip(sol.dqsol.t, sol.dqsol.u)
]

plot(sol.dqsol.t, first.(Δ))
plot!(sol.dqsol.t, last.(Δ))
vline!([t_step,], label="step", c=:black)
vline!([t_fault, t_fault + t_duration], label="fault", c=:black)

# currents
node = 3
plot(
    t_array,
    sol.(t_array, node, :iabs),
    label="i_$node"
)
vline!([t_step,], label="step", c=:black)
vline!([t_fault, t_fault + t_duration], label="fault", c=:black)
# current set point
I_r = [
    complex(sol(t, node, :P_g), -sol(t, node, :Q_g)) / sol(t, node, :v)
    for t in t_array
]
plot!(t_array, abs.(I_r), label="I_r_$node")


P_err = sol.(t_array, 2:3, :P_err)
Q_err = sol.(t_array, 2:3, :Q_err)

plot(first.(P_err), first.(Q_err))
plot!(last.(P_err), last.(Q_err))
