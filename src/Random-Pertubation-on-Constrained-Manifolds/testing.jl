include("../Random-Pertubation-on-Constrained-Manifolds/RandPertOnConstraindManifolds.jl")
include("../../example_cases/BasinStabilityForPD/pd_basin_stability.jl")
include("../../example_cases/BasinStabilityForPD/plotting.jl")

###############################
#            SET UP           #
###############################

op = find_steady_state(powergrid)
rpg = rhs(powergrid)

P, Q = SimplePowerFlow(powergrid, op)

ω_idx = findall(:ω .∈ symbolsof.(powergrid.nodes))
θ_idx = findall(:θ .∈ symbolsof.(powergrid.nodes))

throws = 1000
node = 4

g(θ,i,u) = Q[node] - 3 * i * exp(-2im * θ) * u # Constraint, Power Flow Analysis
z = [op[node,:θ], op[node,:iabs],op[node,:v]]
symbol_list = [:θ, :i, :u]

##################################################

dz = RandPertWithConstrains(g, symbol_list, z, op, node,[-π, π],[-100,100], Uniform)

PerturbedState = copy(op)
PerturbedState[node, :θ] = dz[1]
PerturbedState[node, :v] = dz[3]

#################################
#         PLOTTING              #
#################################

result = sim(powergrid, PerturbedState, (0.0, 20.0), true, rpg)

plot_res(result, powergrid, 3)

plot(PowerGridSolution(result.dqsol, powergrid), θ_idx, :θ, legend = false)

plot_distribution(g, symbol_list, z, op, Uniform, "Optim", node, throws)

plot_distance(g, symbol_list, z, op, Uniform, "Optim", node, throws)
