using LightGraphs
using Test
using Revise,ABMEv
using UnPack
myspace = (GraphSpace(SimpleGraph(10,10)),RealSpace{1,Float64}())
myagents = [Agent(myspace,ancestors=true,rates=true) for i in 1:10]
d(X,Y) = gaussian(X[1],Y[1],1)
b(X,Y) = gaussian(X[1],0,1)
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
@test size(s.agentarray,2) == 2

#TODO: try with callbacks
