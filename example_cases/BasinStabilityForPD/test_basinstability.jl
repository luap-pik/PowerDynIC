include("pd_basin_stability.jl")

line_list = []
append!(line_list, [StaticLine(from = 1, to = 2, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 2, to = 3, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 3, to = 4, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 4, to = 5, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 1, to = 5, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 5, to = 6, Y = -1im / 0.02)])
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


node_list = []
append!(node_list, [PQAlgebraic(S = 0.9, Y_n = 1)])
append!(node_list, [SlackAlgebraic(U = 1, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = -0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_d_dash = 0.111, X_d = 0.1)])
append!(node_list, [PQAlgebraic(S = -0.5, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = -0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_d_dash = 0.111, X_d = 0.1)])
append!(node_list, [PQAlgebraic(S = 0.9, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = 0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_d_dash = 0.111, X_d = 0.1)])
append!(node_list, [PQAlgebraic(S = 0.9, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = 0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_d_dash = 0.111, X_d = 0.1)])
append!(node_list, [PQAlgebraic(S = 0.9, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 3.318, P = 0.6337, D = 0.1, Ω = 50, E_f = 0.5, T_d_dash = 8.690, X_d_dash = 0.111, X_d = 0.1)])

powergrid = PowerGrid(node_list,line_list)
operationpoint = find_operationpoint(powergrid)

node = 7
endtime = 20
guess1 = 29.999995
guess2 = 30.0

result1 = simulate(
    Perturbation(node, :ω, Inc(guess1)),
    powergrid,
    operationpoint,
    timespan = (0.0, endtime),
)

result2 = simulate(
    Perturbation(node, :ω, Inc(guess2)),
    powergrid,
    operationpoint,
    timespan = (0.0, endtime),
)
