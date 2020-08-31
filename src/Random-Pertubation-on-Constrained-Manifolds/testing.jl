include("../Random-Pertubation-on-Constrained-Manifolds/RandPertOnConstraindManifolds.jl")
include("../../example_cases/BasinStabilityForPD/pd_basin_stability.jl")
include("system.jl")

###############################
#            SET UP           #
###############################

op = find_operationpoint(powergrid; sol_method = :rootfind)
pg = powergrid

g_testing(u, resid) = constraint_equations(rhs(pg))

cb = ManifoldProjection(g_testing);
prob = ODEProblem(rhs(pg), x,(0.0, 20.0));

sol = solve(prob, Rosenbrock23(), save = true,callback = cb)
#state = PowerGridSolution(sol, pg)
Plots.plot(sol)

x = RandomWalkManifold(rhs(pg), op.vec, [-1.0,1.0], Uniform, nsteps = 1)
PerturbedState = State(pg, x)
PerturbedState[1,:v]


result = sim(pg, PerturbedState, (0.0,20.0), false, rhs(pg))
plot_res(result, pg , 2)

x = AmbientForcing(rhs(pg), op.vec,[-1.5, 1.5], Uniform, 5, 2)
PerturbedState = State(pg, x)

for j in 1:6
    for i in 1:5
        x = AmbientForcing(rhs(pg), op.vec,[-1.5, 1.5], Uniform, 5, j)
        PerturbedState = State(pg, x)

        result = sim(pg, PerturbedState, (0.0,20.0), true, rhs(pg))
        plot_res(result, pg, 1)
    end
end
