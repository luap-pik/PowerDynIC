include("pd_basin_stability.jl");
include("system.jl");

bs, var = BasinStability(powergrid, 28, 10, [-2,2], Uniform, tau_max = 3);

op = find_operationpoint(powergrid)
sol = sim(powergrid, op, (0.0, 20), true, rhs(powergrid))

gplt = plot_powergrid(powergrid);
array = RemoveNodes(rhs(powergrid), powergrid);

#Plotting
PlotBasinStability(array, bs, var, labtext = "Ambient Forcing");
png("BasinStability_Ambient");
