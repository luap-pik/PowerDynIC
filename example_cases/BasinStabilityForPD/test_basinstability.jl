include("../../src/Random-Pertubation-on-Constrained-Manifolds/system.jl")
include("pd_basin_stability.jl")

bs, var = BasinStability(powergrid,28,:ω,10, [-1000.0,1000.0])
