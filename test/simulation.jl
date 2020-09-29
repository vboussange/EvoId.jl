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
w = World(myagents,myspace,p,0.)

s = Simulation(w)
