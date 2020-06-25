include("../Random-Pertubation-on-Constrained-Manifolds/RandPertOnConstraindManifolds.jl")
include("../../example_cases/BasinStabilityForPD/pd_basin_stability.jl")
op = find_steady_state(powergrid)
rpg = rhs(powergrid)

P, Q = SimplePowerFlow(powergrid, op)

ω_idx = findall(:ω .∈ symbolsof.(powergrid.nodes))
θ_idx = findall(:θ .∈ symbolsof.(powergrid.nodes))

PerturbedState, dz = RandPertWithConstrains(powergrid, op, 5,[-π, π],[-100,100] ,Q[3])

result = sim(powergrid, PerturbedState, (0.0, 20.0), true, rpg)

plot_res(result, powergrid, 3)

plot(PowerGridSolution(result.dqsol, powergrid), θ_idx, :θ, legend = false)

#ND.jl gives us an ODEFunction so we have to construct one manually
#rpg_ode = rhs(powergrid)

##
throws = 1000
node = 4
v_vec = Array{Float64}(undef, throws)
θ_vec = Array{Float64}(undef, throws)
ω_vec = Array{Float64}(undef, throws)
i_vec = Array{Float64}(undef, throws)

for i in 1:throws
        PerturbedState, dz  =  RandPertWithConstrains(powergrid, op, node,[-π, π],[-100,100] ,Q[node])
        v_vec[i] = abs(dz[3])
        θ_vec[i] = abs(dz[1])
        ω_vec[i] = PerturbedState[node,:ω]
        i_vec[i] = PerturbedState[node,:iabs]
end

scatter(θ_vec, v_vec, marker = (:circle, 8), label = "Optim Projection")
xlabel!("θ [rad]")
ylabel!("v [p.u.]")

png("Pertubations_Optim_neu_4")



op[3,:φ]
a = powergrid.nodes[2]

println(a.H)
#=
rr = RootRhs(rhs(powergrid))
rr = (dx,x) -> rpg(dx,x, nothing, 0)
xtol = 1E-9
nl_res = nlsolve(rr, op.vec)
dx = similar(nl_res.zero)

# Define a differential equation system
@parameters t
@variables θ(t) i(t) u(t) g(t)
@derivatives D'~t

# Lorenz attractor
eqs = [g ~ Q[1] - 3 * i * u * exp(-2im * θ)]
de = ODESystem(eqs)

Jg = transpose(calculate_jacobian(de))

println(Jg)
=#
