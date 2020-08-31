using PowerDynamics
using ForwardDiff
using Distributions
using DifferentialEquations
using LinearAlgebra
using Optim

"""
A method to find random pertubations in the case of systems with algebraic
constraints. There is then a set of valid initial conditions which lies on a so
called constrained manifold.
The system must still lie on this manifold after the random pertubation.
This could be applied to Basin Stability Methods in the case of complexer
Generator Methods which contain algebraic constraints.
"""
function RandomWalkManifold(f::ODEFunction, z,
                                dist_args, dist::UnionAll,
                                nsteps)
        g = constraint_equations(f)
        for i in 1:nsteps
            Jacg = ForwardDiff.jacobian(g, z)

            kerJacg = nullspace(Jacg) # kernel of the Jacobian, returns orthogonal vectors

            Frand = rand(dist(dist_args...), size(kerJacg,2))

            dz = kerJacg * Frand

            z = ProjOptim(dz, g) # Calculates the projection onto M using Optim
        end
    return z
end

"""
Calculates and returns the constraint equations from an ODEFunction used in
DifferentialEquations.jl.
The system must be in Mass Matrix form meaning: M ẋ = f(x). The constraints can
therefore be easily extracted by finding the diagonal entries of M which are 0.
Inputs:
       f: An ODEFunction in mass matrix form
Outputs:
        g(x): A function which holds all the constraints of f.

"""
function constraint_equations(f::ODEFunction)
    M = f.mass_matrix
    len_M = size(M, 1)
    g_idx = findall([M[i,i] for i in 1:len_M] .== 0)
    g(x) = (dx = similar(x);
            f(dx, x, nothing, 0.0);
            dx[g_idx])
end

"""
"""
function AmbientForcing(f::ODEFunction, z,
          dist_args, dist::UnionAll, tau_max)
    g = constraint_equations(f)

    M = f.mass_matrix
    len_M = size(M, 1)

    Frand = rand(dist(dist_args...), len_M)
    tau_rand = rand(Uniform(0, tau_max))

    prob = ODEProblem(ambient_forcing_ODE, z, (0.0, tau_rand), (g, Frand))
    sol = solve(prob)
    return sol[end]
end

"""
Perturbed only the variable associated to the node
"""
function AmbientForcing(f::ODEFunction, z,
          dist_args, dist::UnionAll, tau_max, node::Int64)
    g = constraint_equations(f)
    node_vars = NodeToVarIdx(f, node)

    M = f.mass_matrix
    len_M = size(M, 1)

    Frand = zeros(len_M)
    Frand[node_vars] = rand(dist(dist_args...), length(node_vars))
    tau_rand = rand(Uniform(0, tau_max))

    prob = ODEProblem(ambient_forcing_ODE, z, (0.0, tau_rand), (g, Frand))
    sol = solve(prob)
    return sol[end]
end

"""
Takes a ODEFunction from PowerDynamics and returns the indexes of the dynamical
variables associated to a node.
Inputs:
        f: The rigthhand side of a PowerGrid object
        node: Variable indexes of this node are found
Outputs:
        vars_array: Array containg the indexes of the variables of node

"""
function NodeToVarIdx(f::ODEFunction, node::Int64)
    M = f.mass_matrix
    len_M = size(M, 1)
    var_array = []
    for i in 1:len_M
        str = string(f.syms[i])
        idx = findlast('_', str)

        node_num = parse(Int64, str[idx + 1])

        if node_num > node
            return var_array
        end
        if node_num == node
            append!(var_array, i)
        end
    end
    return var_array
end

"""
Calculates the orthogonal projection on a subspace N. The basis of N does not
have to be a orthonormal basis. The matrix (A^T*A)^-1 recovers the norm.
Inputs:
      A: Matrix containg the basis vectors of a subspace N as columns
"""
Proj_N(A) = A * inv(transpose(A) * A) * transpose(A)

"""
Takes a constraint function g and a random value Frand form the ambient and
calculates the projection on to the tangetial space. This gives a manifold
preserving version of any dynamic.
    ż = Pn(ker(Jg)) * Frand
This is an ODE which will be solved using DifferentialEquations.jl
Inputs:
      p[1] = g: The function defining the manifold
      p[2] = Frand: A random constant vector from ambient space
      u0: Inital condition
      t: Integration time
"""
function ambient_forcing_ODE(u, p, t)
    g, Frand = p
    Jacg = ForwardDiff.jacobian(g, u)
    N = nullspace(Jacg)
    du = Proj_N(N) * Frand
end


"""
Calculates the projection back onto the manifold which is given by g.
The function p(z) = b * ∑|g(z)| + a * |dz - z|  is minimized in order to find
the closest value to dz which still fullfills the constraints.
The minimum of p(z) is found by using Optim.
Inputs:
        dz: A value from the tangetial space at some point in the manifold
        g: The function defining the manifold
        a,b: Weigthing Factors
Outputs:
        minimum: Minimum of the optimization process by Optim
"""
function ProjOptim(dz, g::Function; a = 1.0, b = 10.0)
    Proj(z) = b * sum(abs.(g(z))) +
              a * norm(dz - z)

    proj = optimize(Proj, dz, BFGS())

    minimum = proj.minimizer

    return minimum
end
