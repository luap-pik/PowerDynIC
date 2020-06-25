include("../../src/Random-Pertubation-on-Constrained-Manifolds/system.jl")
include("pd_basin_stability.jl")

#bsode, varode = BasinStability(powergrid,28,:ω,100, [-π,π], [-100.0,100.0], ProjectionMethod = "ODE")
#bsoptim, varoptim = BasinStability(powergrid,28,:ω,100, [-π,π], [-100.0,100.0], ProjectionMethod = "Optim")

bs, var = BasinStability(powergrid,28.0,:ω,1, [-π,π], [-100.0,100.0], ProjectionMethod = "Optim", dae = true)
odebs, varode = BasinStability(powergrid,28,:ω,1, [-π,π], [-100.0,100.0], ProjectionMethod = "Optim")


array = RemoveNodes(powergrid,:θ)
PlotBasinStability(array, bs, var, labtext = "Optim Projection, ODE")
PlotBasinStability!(array, odebs, varode, labtext = "Optim Projection, DAE")
png("BasinStability_Plot_ode_dae")


using DelimitedFiles
writedlm("var_ode.txt", varode)
writedlm("var_optim.txt", varoptim)


using GraphPlot

using LightGraphs
using Cairo, Compose
draw(PNG("powergrid.png", 16cm, 16cm), gplot(powergrid.graph))

vec = powergrid.nodes |> collect

types = typeof.(vec)
