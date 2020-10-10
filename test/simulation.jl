using LightGraphs
using Test
using Revise,ABMEv
using UnPack
myspace = (GraphSpace(SimpleGraph(10,10)),RealSpace{1,Float64}())
myagents = [Agent(myspace,ancestors=true,rates=true) for i in 1:10]
d(X,Y,t) = gaussian(X[1],Y[1],1)
b(X,Y,t) = gaussian(X[1],0,1)
D = (Int64(1),Float64(1.))
mu = [1,1]
NMax = 100
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax
w0 = World(myagents,myspace,p,0.)

# Basic test
s = Simulation(w0)
typeof(s)
@test get_tend(s) ≈ 0.
@test get_size(s) ≈ 1

# adding an entry to sim
newa = give_birth(1,w0)
addAgent!(w0,newa)
@test isnothing(add_entry!(s,w0))
@test get_size(s) == 2

#TODO: try with callbacks
###################################
#######CALLBACKS###################
###################################

myspace = (RealSpace{1,Float64}(),)
sigma_K = .9;
sigma_a = .7;
K0 = 1000;
b(X) = gaussian(X[1],0.,sigma_K)
d(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
D = (1e-2,)
mu = [.1]
NMax = 10000
tend = 1.5
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax

myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
w1 = copy(w0)
cb = (names = ["gamma_div"], agg = Function[w -> var(Float64.(get_x(w,1)))])
eltype(cb.agg)
@time sim = run!(w1,Gillepsie(),tend,cb=cb,dt_saving = .1)

@test typeof(sim["gamma_div"]) <: Vector

@test typeof(get_world(sim,get_size(sim))) <: World
@test typeof(sim[2]) <: Vector

## testing plot
using Plots
plot(get_tspan(sim),sim["gamma_div"])
