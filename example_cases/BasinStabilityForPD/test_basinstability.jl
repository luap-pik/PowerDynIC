include("../../src/Random-Pertubation-on-Constrained-Manifolds/system.jl")
include("pd_basin_stability.jl")

BasinStability(powergrid,20,:ω,20, [-25.0,25.0])
