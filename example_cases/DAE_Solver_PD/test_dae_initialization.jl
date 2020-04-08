include("dae_initialization.jl")
include("ThirdOrderEq.jl")

line_list = []
append!(line_list, [StaticLine(from = 1, to = 2, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 2, to = 3, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 1, to = 3, Y = -1im / 0.02)])

node_list = []
append!(node_list, [SlackAlgebraic(U = 1, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = -0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_d_dash = 0.111, X_d = 0.1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = -0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_d_dash = 0.111, X_d = 0.1)])

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
