include("../../src/Random-Pertubation-on-Constrained-Manifolds/system.jl")
include("pd_basin_stability.jl")

bs, var = BasinStability(powergrid,28,:ω,100, [-10.0,10.0], [-π,π])

using GraphPlot
using LightGraphs
using Cairo, Compose
draw(PNG("powergrid.png", 16cm, 16cm), gplot(powergrid.graph))

vec = powergrid.nodes |> collect

types = typeof.(vec)
