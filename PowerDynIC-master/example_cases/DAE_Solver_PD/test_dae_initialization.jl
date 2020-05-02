include("dae_initialization.jl")
include("../../components/ThirdOrderEq.jl")

include("system.jl")

powergrid = PowerGrid(node_list,line_list)
operationpoint = find_operationpoint(powergrid)

timespan = (0.0, 30)

# I choose du0 = 0 because I assume that we start from a stable fixpoint for now
# but this should definitely change later on....
# Seems like a valid DAEProblem is found. The function is of type nd_ODE_ODE{...}
DAEProblem(dae_rhs(powergrid),zeros(systemsize(powergrid)),operationpoint.vec,timespan)

# The function looks the same as in the DAEProblem
ODEProblem{iipfunc}(rhs(powergrid),operationpoint.vec,timespan)

solve(powergrid, operationpoint,timespan)

# Fails because there is no method for nd_ODE_ODE{...}.
# But exactly this type seems to work in the case of ODE solvers...
dae_solve(powergrid, operationpoint, timespan)
