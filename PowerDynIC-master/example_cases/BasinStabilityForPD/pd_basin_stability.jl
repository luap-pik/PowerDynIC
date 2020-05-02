using PowerDynamics
using Distributions
using NetworkDynamics
using LaTeXStrings
include("plotting.jl")
include("ThirdOrderEq.jl")


"""
Evaluates the result of a simulation after a pertubation.
Checks if conditions for a stable synchronous state where met.
Conditions implemented so far:
1. Is ω returning to zero? (Realized through checking the deviations from 0)
#2. Did one of the machines get too far away from a stable state? -> Survivability
3. How close are the phases θ to each other? -> Braucht man das?
Inputs:
       pg: Power grid, a graph containg nodes and lines
       result: Solution of a power grid after a pertubation
       nodesarray: Array containg the index of all nodes that are not Slacks
       endtime: Endtime of the simulation
Outputs:
       Returns 1 if the system fullfills the condition for a stable state
       Returns 0 if the condition was not met
"""
function StableState(result::PowerGridSolution,nodesarray,variable,endtime, interval)
    summation = 0
    for node in nodesarray
        #=
        maxdeviation = maximum(result(0:endtime,node,:ω))
        if abs(maxdeviation) > interval[2]               # this could be used for Survivability!
            return 0
        end
        =#
        for time in endtime-10:0.1:endtime
            summation =+ abs.(result(time,node,:ω))
        end
        #summation =+ sum(abs.(result(endtime-10:0.1:endtime,node,:ω))) # this does not seem to work anymore since i updated Pd...
    end
    if summation > 100 * length(nodesarray)
        return 0
    else
        return 1
    end
end

"""
Generates a random pertubation applied to a single node.
It is then checked if a valid initital conditions can be found on the constraind
manifold using the PowerDynamics function find_valid_initial_condition.
This is repeated until valid conditions can be found.
Inputs:
      pg: Power grid, a graph containg nodes and lines
      node: Node index for the perturbed node
      variable:  Internal variable of the node types to be perturbed
      interval: Lower and Upper Bound for the pertubation
outputs:
        initial_state: Initial state for the simulation
        randpertubation: Random number inside interval
"""
function ConstrainedManifold(pg::PowerGrid, node, variable, interval)
    operationpoint = find_operationpoint(pg)
    perturbedstate = operationpoint

    while true
        randpertubation = rand(Uniform(interval[1],interval[2]))
        println(randpertubation)
        perturbedstate[node, variable] = randpertubation
        try
            state = find_valid_initial_condition(pg,perturbedstate.vec)
            return state, randpertubation
        catch
        end
    end
end

"""
Removes all SlackBuses from the calculation  procedure. Returns a list
containg the indexes of all nodes that are not of type SlackAlgebraic.
Inputs:
    pg: Power grid, a graph containg nodes and lines
Outputs:
    nodesarray: Array containg the index of all nodes that are not Slacks
"""
function RemoveSlackBuses(pg::PowerGrid)
    nodesarray = []
    for node = 1:size(pg.nodes)[1]
        if typeof(pg.nodes[node]) == SlackAlgebraic
            continue
        end
        if typeof(pg.nodes[node]) == PQAlgebraic
            continue
        end
        append!(nodesarray,node)
    end
    return nodesarray
end


"""
Produces a plot of the calulated basin stability against the nodes index.
Inputs:
    nodesarray: Array containg the index of all nodes that are not Slacks
    basinstability: Array containing the basin stability of each node
    standarddev: Array of the standard deviations of the basin stability
"""
function PlotBasinStability(nodesarray, basinstability, standarddev)
    plot(
        nodesarray,
        basinstability,
        #yerror = standarddev,
        #title = "...",
        xlabel = "Node Index",
        ylabel = "Basin Stability",
        marker = (:circle, 8, "red"),
        line = (:path, 2,"gray"),
        legend = false,
        grid = false,
        show = true
    )
end

"""
Calculates the basin stability of the nodes in a power grid using
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
The behaviour of the power grid after the pertubation are then calculated.
In the function StableState the number of times the system returns
to the stable synchronous state is counted. If numtries -> ∞:
Then BasinStability ≈ stablecounter / numtries.
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
function BasinStability(pg::PowerGrid, endtime::Int, variable::Symbol,
                        numtries::Int, interval::Vector{Float64})
    @assert interval[2] >= interval[1] "Upper bound should be bigger than lower bound!"

    operationpoint = find_operationpoint(pg)
    nodesarray = RemoveSlackBuses(pg)
    basinstability = []
    standarddev =  []

    for node in nodesarray
        stablecounter = Array{Int64}(undef,numtries)
        for tries = 1:numtries
            println("NODE:",node," TRY:",tries)
            state, pertubation = ConstrainedManifold(pg,node,variable, interval)
            try
                result = simulate(
                    Perturbation(node, variable, Inc(pertubation)),
                    pg,
                    state,
                    timespan = (0.0, endtime),
                )
                stablecounter[tries] = StableState(result,nodesarray,variable,endtime,interval)

                stablecounter[tries] = StableState(result,nodesarray,variable,endtime,interval)
            catch err
                if isa(err, GridSolutionError)
                    stablecounter[tries] = 0 # counted as a unstable solution
                end
            end
        end
        append!(basinstability, mean(stablecounter))
        append!(standarddev, std(stablecounter)) # is this even a useful measure here?
    end

    display(PlotBasinStability(nodesarray, basinstability, standarddev))
    png("BasinStability_Plot")
    return basinstability, standarddev
end
