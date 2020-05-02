"""
Allows the user to set a starting point for a grid in a PowerDynamics simulation.
Inputs:
        pd: Powergrid to be simulated
        v_vec: Vector of voltages at the nodes
        ω_vec: Vector of frequencies at the nodes
        θ_vec: Vector of phases at the nodes
Outputs:
        Startingpoint: Customized State of the grid
"""
using PowerDynamics: PowerGrid, State

function my_startingpoint(pd::PowerGrid,v_vec, ω_vec, θ_vec)
    numnodes = size(pd.nodes)[1]
    @assert size(v_vec)[1] == numnodes "Size of v_vec not the same as number of nodes."
    @assert size(ω_vec)[1] == numnodes "Size of ω_vec not the same as number of nodes."
    @assert size(θ_vec)[1] == numnodes "Size of θ_vec not the same as number of nodes."

    println("Warning! Manually setting the startingpoint might not result
    in a stable fixpoint.")

    startingpoint = State(pd, ones(systemsize(pd))) # just put in a bunch of ones for now
    for i in 1:numnodes
        startingpoint[i,:ω] = ω_vec[i]
        startingpoint[i,:v] = v_vec[i]
        startingpoint[i,:θ] = θ_vec[i]
    end
    return startingpoint
end
