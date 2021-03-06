include("../../example_cases/new node types/plotting.jl")
include("../../example_cases/new node types/ThirdOrderEq.jl")
include("../../example_cases/new node types/ThirdOrderEqGovernorIEEEG1.jl")
include("../../example_cases/new node types/ThirdOrderEqExciterIEEEDC1A.jl")

line_list = []
append!(line_list, [StaticLine(from = 1, to = 2, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 2, to = 3, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 3, to = 4, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 4, to = 5, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 5, to = 6, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 1, to = 6, Y = -1im / 0.02)])

#=
append!(line_list, [StaticLine(from = 6, to = 7, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 7, to = 8, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 8, to = 9, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 9, to = 10, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 10, to = 11, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 1, to = 11, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 2, to = 7, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 1, to = 4, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 2, to = 11, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 2, to = 9, Y = -1im / 0.02)])
=#

node_list = []

append!(node_list, [SlackAlgebraic(U = 1, Y_n = 0)])
#append!(node_list, [PQAlgebraic(S = 0.9, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = -0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_q_dash = 0.103, X_d_dash = 0.111, X_d = 0.1)])
#append!(node_list, [PQAlgebraic(S = -0.5, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = -0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_q_dash = 0.103, X_d_dash = 0.111, X_d = 0.1)])
#append!(node_list, [PQAlgebraic(S = 0.9, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = 0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_q_dash = 0.103, X_d_dash = 0.111, X_d = 0.1)])
#append!(node_list, [PQAlgebraic(S = 0.9, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = 0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_q_dash = 0.103, X_d_dash = 0.111, X_d = 0.1)])
#append!(node_list, [PQAlgebraic(S = 0.9, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = 0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_q_dash = 0.103, X_d_dash = 0.111, X_d = 0.1)])



#=
append!(node_list, [SlackAlgebraic(U = 1, Y_n = 1)])
append!(node_list, [ThirdOrderEqGovernorIEEEGone(H = 5.0, D = 0.1, Ω = 50, E_f = 0.2, T_d_dash = 0.1, X_d_dash = 0.5,  X_q_dash = 0.103, X_d = 0.5, P0 = 0.5, Pmax = 1.03, Pmin = 0, Pup = 1.38, Pdown = -1.25, T1 = 0.82, T2 = 0, T3 = 0.97, K = 20.0)])
append!(node_list, [ThirdOrderEqGovernorIEEEGone(H = 5.0, D = 0.1, Ω = 50, E_f = 0.2, T_d_dash = 0.1, X_d_dash = 0.5,  X_q_dash = 0.103, X_d = 0.5, P0 = 0.5, Pmax = 1.03, Pmin = 0, Pup = 1.38, Pdown = -1.25, T1 = 0.82, T2 = 0, T3 = 0.97, K = 20.0)])
append!(node_list, [ThirdOrderEqGovernorIEEEGone(H = 5.0, D = 0.1, Ω = 50, E_f = 0.2, T_d_dash = 0.1, X_d_dash = 0.5,  X_q_dash = 0.103, X_d = 0.5, P0 = 0.5, Pmax = 1.03, Pmin = 0, Pup = 1.38, Pdown = -1.25, T1 = 0.82, T2 = 0, T3 = 0.97, K = 20.0)])
append!(node_list, [ThirdOrderEqGovernorIEEEGone(H = 5.0, D = 0.1, Ω = 50, E_f = 0.2, T_d_dash = 0.1, X_d_dash = 0.5,  X_q_dash = 0.103, X_d = 0.5, P0 = 0.5, Pmax = 1.03, Pmin = 0, Pup = 1.38, Pdown = -1.25, T1 = 0.82, T2 = 0, T3 = 0.97, K = 20.0)])
append!(node_list, [ThirdOrderEqGovernorIEEEGone(H = 5.0, D = 0.1, Ω = 50, E_f = 0.2, T_d_dash = 0.1, X_d_dash = 0.5,  X_q_dash = 0.103, X_d = 0.5, P0 = 0.5, Pmax = 1.03, Pmin = 0, Pup = 1.38, Pdown = -1.25, T1 = 0.82, T2 = 0, T3 = 0.97, K = 20.0)])

=#


powergrid = PowerGrid(node_list, line_list)

return powergrid
