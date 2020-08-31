using PowerDynamics: simulate, Perturbation, Inc, PowerGrid
using Distributions: Uniform
using Sundials
#using Statistics
using NetworkDynamics
using LaTeXStrings
#using FileIO, JLD2
#using DelimitedFiles

include("../new node types/ThirdOrderEq.jl")
include("../../src/Random-Pertubation-on-Constrained-Manifolds/RandPertOnConstraindManifolds.jl")
#include("../../src/PDpatches.jl")
include("plotting.jl")

"""
Evaluates the result of a simulation after a pertubation.
Checks if conditions for a stable synchronous state where met.
Conditions implemented so far:
1. Is ω returning to zero? (Realized through checking the deviations from 0)
#2. Did one of the machines get too far away from a stable state? -> Survivability
Inputs:
       pg: Power grid, a graph containg nodes and lines
       result: Solution of a power grid after a pertubation
       nodesarray: Array containg the index of all nodes that contain varible
       endtime: Endtime of the simulation
Outputs:
       Returns 1 if the system fullfills the condition for a stable state
       Returns 0 if the condition was not met
"""
function StableState(pg, result::PowerGridSolution,nodesarray,endtime)
    ω_idx = findall(:ω.∈ symbolsof.(pg.nodes))
    for node in  1:length(pg.nodes) # runs over all nodes
        summation = 0
        for time in endtime-10:0.1:endtime
            if node in ω_idx  # checks only nodes which have ω as a symbol
                summation =+ abs.(result(time,node,:ω))
            end
            if summation > 100  # checks if ω returns to 0
                                #this is just an arbitrary value that i have choosen so far
                return 0
            end
            if isapprox(result(time,node, :v) , 0, atol=1e-2) # checks for a short cirtuit
                return 0
            end
        end
    end
    return 1
end

function sim(pg::PowerGrid, x0::State, timespan, force_dtmin::Bool, rpg)
    problem = ODEProblem(rpg,x0.vec,timespan)
    # sets the kwarg force_dtmin to avoid the dt < dtmin error
    if force_dtmin == true
        solution = solve(problem, Rodas5(autodiff=true), force_dtmin = true)
    else
        solution = solve(problem, Rodas5(autodiff=false))
    end
    PowerGridSolution(solution, pg)
end

function dae_sim(rpg_ode, tspan, x0::State,pg)
        function dae_form(res, du, u, p, t)
            du_temp = similar(du)
            rpg_ode.f(du_temp, u, p, t)
            @. res = du - du_temp
        end
        u0 = x0.vec
        du0 = similar(u0)
        p = nothing
        diff_vars = (diag(rpg_ode.mass_matrix) .== 1)
        rpg_dae = DAEFunction{true, true}(dae_form, syms = rpg_ode.syms)
        dae = DAEProblem(rpg_dae, du0, u0, tspan, p; differential_vars=diff_vars)
        sol = solve(dae, DABDF2(autodiff=false), initializealg=BrownFullBasicInit(), force_dtmin = true, matiters= 1e7,save_everystep=false)
        #sol = solve(dae, IDA(linear_solver=:GMRES))
        return PowerGridSolution(sol,pg)
end


"""
Removes all nodes from the from the calculation procedure that do not contain
variable as a symbol.
Returns a list containg the indexes of all the remaning nodes.
Inputs:
    pg: Power grid, a graph containg nodes and lines
Outputs:
    nodesarray: Array containg the index of all nodes with the symbol varibale
"""
function RemoveNodes(rpg, pg::PowerGrid)
    constraint_idx = findall(0 .== [rpg.mass_matrix[i, i] for i in 1:systemsize(pg)])
    nodesarray = []
    for i in constraint_idx
        symbol_str = string(rpg.syms[i])
        num_idx = findlast(isequal('_'), symbol_str) + 1
        node_num = parse(Int64, symbol_str[num_idx:end])
        if (node_num in nodesarray) == false
            append!(nodesarray, node_num)
        end
    end
    return nodesarray
end


"""
Calculates the single node basin stability of the nodes in a power grid using
PowerDynamics.jl.
The calculations were performed on the basis of the three following papers:
1. Menck, P., Heitzig, J., Marwan, N. et al. How basin stability complements
   the linear-stability paradigm. Nature Phys 9, 89–92 (2013)
2. Menck, P., Heitzig, J., Kurths, J. et al. How dead ends undermine
   power grid stability. Nat Commun 5, 3969 (2014).
3. C Mitra, A Choudhary, S Sinha, J Kurths, and R V Donnerk,
   Multiple-node basin stability in complex dynamical networks
   Physical Review E, (2017)

To calculate the basin stability "numtries" random pertubations are
drawn according to the interval.
The behaviour of the power grid after the pertubation is then calculated.
In the function StableState the number of times the system returns
to the stable synchronous state is counted.
If numtries -> ∞: Then BasinStability ≈ stablecounter / numtries.
Inputs:
    pg: Power grid, a graph containg nodes and lines
    endtime: Endtime of the simulation
    variable: Internal variable of the node types to be perturbed
    numtries: Number of random pertubations within interval to be evaluated
    interval: Lower and Upper Bound for the pertubation
Outputs:
    basinstability: Array conating the basinstability of each node
    standarddev: Array of the standard deviations of the basin stability
"""
function BasinStability(pg::PowerGrid, endtime,
                        numtries::Int, dist_args, dist::UnionAll;
                        tau_max)
    num_errors = 0
    op = find_operationpoint(pg; sol_method = :rootfind)
    rpg = rhs(pg)

    nodesarray = RemoveNodes(rpg, pg)
    basinstability = []
    standarddev =  []

    #state_file = File(format"JLD2","State_File.jld2")
    #save(state_file, "operationpoint", op)
    for node in nodesarray
        stablecounter = Array{Int64}(undef,numtries)

        for tries = 1:numtries
            println("NODE:",node," TRY:",tries)
            z_new = AmbientForcing(rpg, op.vec, dist_args, dist, tau_max, node)

            PerturbedState = State(pg, z_new)

            #column_name = string(node) * string(tries)
            #save(state_file, column_name, PerturbedState)
            λ, stable = check_eigenvalues(powergrid, PerturbedState)
            #=
            if stable != true
                println("Positive eigenvalue detected. The system might be unstable.")
                try
                    result = sim(pg, PerturbedState, (0.0, endtime),true, rpg)
                    stablecounter[tries] = StableState(result, nodesarray, variable, endtime, interval)
                    if stablecounter[tries] == 0
                        display(plot_res(result, pg, node))
                        fn = "node$(node)_try$(tries)"
                        png(fn)
                    end
                catch err
                    if isa(err, GridSolutionError)
                        stablecounter[tries] = 0 # counted as an unstable solution
                        num_errors += 1
                    end
                end
            else
            =#
            result = sim(pg, PerturbedState, (0.0, endtime), true, rpg)
            stablecounter[tries] = StableState(pg, result, nodesarray, endtime)
            #if stablecounter[tries] == 0
            display(plot_res(result, pg, node))
            fn = "node$(node)_try$(tries)"
            png("/home/anna/Documents/MASTER_ARBEIT/random stuff/"*fn)
        end
        append!(basinstability, mean(stablecounter))
        append!(standarddev, std(stablecounter)) # is this even a useful measure here?
    end
    println("Percentage of failed Simualtions:  ",100 * num_errors/(numtries * length(nodesarray)))
    display(PlotBasinStability(nodesarray, basinstability, standarddev))
    png("BasinStability_Plot")
    #writedlm("bs.txt", basinstability)
    return basinstability, standarddev
end
