using DifferentialEquations
using Plots
using DiffEqBase
using PowerDynamics
using NetworkDynamics
using NetworkDynamics: collect_ve_info

const iipfunc = true # is in-place function

"""
The almost the same as the original NetworkDynamics-function. Only ODE is switched to DAE
"""
function dae_network_dynamics(vertices, edges, graph; x_prototype = zeros(1)) #where {T <: ODEVertex, U <: ODEEdge}
    v_dims, e_dims, symbols_v, symbols_e, mmv_array, mme_array = collect_ve_info(vertices, edges, graph)

    x_array = similar(x_prototype, sum(v_dims) + sum(e_dims))

    v_array = view(x_array, 1:sum(v_dims))
    e_array = view(x_array, sum(v_dims)+1:sum(v_dims)+sum(e_dims))

    symbols = vcat(symbols_v, symbols_e)

    graph_stucture = GraphStruct(graph, v_dims, e_dims, symbols_v, symbols_e)

    graph_data = GraphData(v_array, e_array, graph_stucture)

    nd! = nd_ODE_ODE(vertices, edges, graph, graph_stucture, graph_data)

    mass_matrix = construct_mass_matrix(mmv_array, graph_stucture)

    #DAEFunction{iipfunc,true}(nd!; mass_matrix = mass_matrix, syms=symbols)

    # mass_matrix is not a supported kwarg here... ill look into that later
    DAEFunction{iipfunc,true}(nd!;syms=symbols)
end

"""
Almost the same as the PowerDynamics function but calls dae_network_dynamics
"""
function dae_rhs(pg::PowerGrid)
    dae_network_dynamics(map(construct_vertex, pg.nodes), map(construct_edge, pg.lines), pg.graph)
end



"""
Almost the same as the PowerDynamics-function, Switched to a DEA Problem!
"""
function dae_solve(pg::PowerGrid, x0, timespan)
    problem = DAEProblem(dae_rhs(pg),zeros(systemsize(pg)),x0.vec,timespan)
    solution = solve(problem,DABDF2(),initializealg = ShampineCollocationInit())
    PowerGridSolution(solution, pg)
end


"""
Calls dae_solve
"""
function dae_simulate(p::Perturbation, powergrid, x0; timespan)
    dae_solve(powergrid, p(x0), timespan);
end
