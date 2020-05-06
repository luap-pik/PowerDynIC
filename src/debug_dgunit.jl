using PowerDynamics
using OrdinaryDiffEq
using DifferentialEquations: DynamicSS
using Sundials
using Plots

const base_power = 1E6 # W
const base_voltage = 20E3 # V
const base_current = base_power / base_voltage # A
const base_admittance = base_power / base_voltage^2 # 1/Ω
const ω = 2 * π * 50.0 # Hz

# some additional functions that might develop into new PD.jl features
include("PDpatches.jl")
# custom node type
include("components/ThirdOrderEq.jl")
include("components/DGUnit.jl")

function testbed(device, U, S, Y)
    busses = [SlackAlgebraic(;U=U), device, PQAlgebraic(;S=S)]
    lines = [StaticLine(;from=1, to=2, Y=Y), StaticLine(;from=2, to=3, Y=Y)]
    return PowerGrid(busses, lines)
end

device = ThirdOrderEq(H = 3.318, P = 0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_q_dash = 0.103, X_d_dash = 0.111, X_d = 0.1)

# no background load
device = DGUnit(
    P_bkgrnd = 0.,
    Q_bkgrnd = 0.,
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
    Vac_ref = 20E3  / base_voltage, #* sqrt(2/3),
    V_dead = 0.1 * 20E3  / base_voltage, #* sqrt(2/3),
    k_FRT = 2.0 / base_current, #/ base_current,
    imax = 1e6/20e3 / base_current / sqrt(3), # MW/kV
    K_pll = 0.1 * base_voltage,
    Y_shunt = t -> 0.#-0.1 * SmoothStep(t, 3.15, 3.; order=100) #0.,
)

pg = testbed(device, 1+0im, -0.6337, -50im)
p = ([(.1, 0., false, false) for _ in 1:length(pg.nodes)], nothing)
state(vec) = State(pg, vec)
solution(vec) = PowerGridSolution(vec, pg)

op_guess = initial_guess(pg)

ode = rhs(pg)

# dae = rhs_dae(pg)
# prob = DAEProblem(dae, similar(op.vec), op.vec, (0., 10.), nothing; differential_vars=differential_vars(ode))
# sol = solve(prob, IDA())
# sol = solve(prob, DABDF2(), initializealg=ShampineCollocationInit())

prob = ODEProblem(ode, op_guess, (0., 10.), p)
op = solve(SteadyStateProblem(prob), DynamicSS(Rodas5(); abstol=1e-8,reltol=1e-6,tspan=Inf)) |> state

# TODO: debug perturb
# TODO: why is there no change???
u0 = copy(op.vec)
v = op[2, :v] * exp(1im * π/3) # phase shift
u0[[3, 4]] .= [real(v), imag(v)]

prob = ODEProblem(ode, op_guess, (0., 10.), p)
sol = solve(prob, Rodas4()) |> solution

plot(sol, :, :v)
