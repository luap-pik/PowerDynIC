using PowerDynamics
using OrdinaryDiffEq

##
dir =  @__DIR__
include("$dir/real_comps.jl")



## parameters from Raphael

omega_0 = 0.
VBase = 320e3
WBase = 1e9
ratio = 10
kappa = atan(ratio)
lengthLine = 100.
res = 0.03 * lengthLine * WBase / VBase^2 # VBase^2/WBase
adm = 1 / (res + im * ratio * res)

## other components

L1 = LineReal(; from = 1, to = 2, Y = adm)
L2 = LineReal(; from = 2, to = 3, Y = adm)
S = SlackReal(; U = complex(1.)) # 1pu
T = ThirdOrderReal(H = 5.0, P = 1.0, D = 0.1, Ω = 50, E_f = 1., X_d_dash = 0.111, T_d_dash = 0.1, X_q_dash = 0.103, X_d = 0.1)
PQ = PQReal(; P = -0.5, Q=-0.1)

##

pg = PowerGrid([S, T, PQ], [L1, L2])

op = find_operationpoint(pg, sol_method=:rootfind)

##

rpg = rhs(pg)

ode = ODEProblem(rpg, op.vec, (0., 10.))

##

using ModelingToolkit

# https://mtk.sciml.ai/stable/tutorials/modelingtoolkitize/

mtk_sys = modelingtoolkitize(ode)

jac = generate_jacobian(mtk_sys)[1]

J = eval(jac)

J(op.vec, nothing, 0.) |> eigvals



##
using ForwardDiff: jacobian

constraint_idx = findall( 0 .== [rpg.mass_matrix[i, i] for i in 1:systemsize(pg)])

f!(dx, x) = rpg(dx, x, nothing, 0.)
f(x) = (dx = similar(x); rpg(dx, x, nothing, 0.); dx)

J = jacobian(f, op.vec)




