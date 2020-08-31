include("../Random-Pertubation-on-Constrained-Manifolds/RandPertOnConstraindManifolds.jl")
include("../../example_cases/BasinStabilityForPD/plotting.jl")

function rober(du,u,p,t)
  y₁,y₂,y₃ = u
  du[1] = -0.04 * y₁ + 1e4 * y₂ * y₃
  du[2] =  0.04 * y₁ - 1e4 * y₂ * y₃ - 3e7 * y₂^2
  du[3] =  y₁ + y₂ + y₃ - 1
  nothing
end

M = [1. 0  0
     0  1. 0
     0  0  0];

ode_rober = ODEFunction(rober,mass_matrix=M);
u0 = [1.0,0.0,0.0]

x = AmbientForcing(ode_rober, u0, [-10,10], Uniform, 5)

#######################################################
#        Drawing Random Inital Conditions             #
###############################################

y1(y2,y3) = -y2 -y3 + 1.0
var_list = [:y1, :y2, :y3]

###################################
#            PLOTTING             #
###################################
scene =  distribution_Makie(ode_rober, y1, var_list, u0, [-100,100], Uniform, "RandomWalk", 5)

scene.center = false
Makie.save("RobertsonRandomWalk.png", scene)

display(scene)
