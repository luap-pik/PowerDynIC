"""
The Robertson Equations in DAE. Used to describe chemical rate equations.
Taken from the DifferentialEquations.jl documentation.
https://diffeq.sciml.ai/stable/tutorials/dae_example/
"""
function Robertson(out,du,u,p,t)
  out[1] = - 0.04u[1]              + 1e4*u[2]*u[3] - du[1]
  out[2] = + 0.04u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
  out[3] = u[1] + u[2] + u[3] - 1.0
end

u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
tspan = (0.0,100000.0)


diff_vars = [true,true,false]
prob = DAEProblem(Robertson,du₀,u₀,tspan,differential_vars = diff_vars)

sol = solve(prob,IDA())

plot(sol, xscale=:log10, tspan=(1e-6, 1e5), layout=(3,1))


gRobertson(y1, y2, y3) = y1 + y2 + y3 - 1.0

y1(y2,y3) = -y2 -y3 - 1.0
node = 1

var_list = [:y1, :y2, :y3]

dz = RandPertWithConstrains(gRobertson, var_list, u₀, op, node,[-10, 10],[-10,10], Uniform)

extrem_y, extrem_z = plot_distribution(gRobertson, var_list, u₀, op, Uniform, "Optim", node)

test = [(y1(0,extrem_y[1]),0,extrem_y[1]), (y1(0,extrem_y[2]),0,extrem_y[2]),(y1(extrem_z[1],0),extrem_z[1],0), (y1(extrem_z[2],0),extrem_z[2],0)]

surface!(test, c = :blues)
