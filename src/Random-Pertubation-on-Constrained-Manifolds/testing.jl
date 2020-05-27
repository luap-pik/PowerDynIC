include("../Random-Pertubation-on-Constrained-Manifolds/RandPertOnConstraindManifolds.jl")
include("../../example_cases/BasinStabilityForPD/pd_basin_stability.jl")
op = find_operationpoint(powergrid)
rpg = rhs(powergrid)
P, Q = SimplePowerFlow(powergrid, op)

ω_idx = findall(:ω .∈ symbolsof.(powergrid.nodes))

θ_idx = findall(:θ .∈ symbolsof.(powergrid.nodes))

PerturbedState, dz = RandPertWithConstrains(powergrid, op, 5, [-100,100], Q[5])

check_eigenvalues(powergrid, PerturbedState)

result = sim(powergrid, PerturbedState, (0.0,100.0),false, rpg)

plot_res(result, powergrid, 3)

plot(PowerGridSolution(result.dqsol, powergrid), θ_idx, :θ, legend = false)

##
throws = 100
node = 5
v_vec = Array{Float64}(undef, throws)
θ_vec = Array{Float64}(undef, throws)
ω_vec = Array{Float64}(undef, throws)
for i in 1:throws
        PerturbedState, dz  = RandPertWithConstrains(powergrid, op, node, [-100,100])
        v_vec[i] = abs(dz[3])
        θ_vec[i] = abs(dz[1])
        ω_vec[i] = PerturbedState[node,:ω]
end

scatter(θ_vec, v_vec, marker = (:circle, 8))
xlabel!("theta [rad]")
ylabel!("v [p.u.]")
