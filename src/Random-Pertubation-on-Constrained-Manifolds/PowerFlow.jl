"""
Caluclates and returns the admittance laplacian of the powergrid pg.
The notation was choosen in accordance to: https://en.wikipedia.org/wiki/Nodal_admittance_matrix
Inputs:
        pg: Powergrid, a graph contaning nodes and lines
"""
function admittance_laplacian(pg::PowerGrid)
    Y =  zeros(Complex, length(pg.nodes), length(pg.nodes))
    y_sum = zeros(Complex, length(pg.nodes))
    for i in 1:length(pg.lines)
        From = pg.lines[i].from
        To = pg.lines[i].to
        Y[From,To] = -1 * pg.lines[i].Y
        y_sum[From] += pg.lines[i].Y
    end
    for j in 1:length(pg.nodes)
        Y[j,j] = pg.nodes[j].Y_n + y_sum[j]
    end
    return Y
end


"""
A Simple Power Flow Calculation for PowerDynamics.
The calculation is only correct for nodes where the phase θ,
voltage magnitude u are know. For powergrids contaning eg. SlackAlgebraic the
calculation can not be performed because the phases θ are unknown.
A Power-flow study would need to be performed before hand.
Inputs:
    pg: PowerGrid, a graph containg nodes and lines
    op: State, Operationpoint (or any other valid point) of the powergrid pg
Outputs:
    P: Vector, Active Power at the nodes
    Q: Vector, Reactive Power at the nodes
"""
function SimplePowerFlow(pg::PowerGrid, op::State)
    #if SlackAlgebraic ∈ pg.nodes .|> typeof
    #    throw(OperationPointError("Found a SlackAlgebraic in the powergrid."))
    #end
    Y = admittance_laplacian(pg)
    G = real(Y)
    B = imag(Y)

    P = zeros(Complex, length(pg.nodes))
    Q = zeros(Complex, length(pg.nodes))

    θ_idx = findall(:θ .∈ symbolsof.(pg.nodes))
    for i in θ_idx
        Vi = op[i, :v] # Voltage Magnitude of node i
        for k in θ_idx
            Vk = op[k,:v] # Volatge Magnitude of node k
            θik = op[i, :θ] - op[k, :θ] # Phase difference between node i and k
            P[i] += Vk * (G[i,k] * sin(θik) + B[i,k] * sin(θik))
            Q[i] += Vk * (G[i,k] * cos(θik) - B[i,k] * cos(θik))
        end
        P[i] *= Vi
        Q[i] *= Vi
    end
    return P, Q
end
