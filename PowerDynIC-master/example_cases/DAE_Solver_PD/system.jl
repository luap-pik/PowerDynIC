line_list = []
append!(line_list, [StaticLine(from = 1, to = 2, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 2, to = 3, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 1, to = 3, Y = -1im / 0.02)])

node_list = []
append!(node_list, [SlackAlgebraic(U = 1, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = -0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_q_dash = 0.103, X_d_dash = 0.111, X_d = 0.1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = -0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_q_dash = 0.103, X_d_dash = 0.111, X_d = 0.1)])
