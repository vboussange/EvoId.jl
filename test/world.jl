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



@test eltype(myagents) <: AbstractAgentM
@test typeof(myagents) <: Vector{A} where {A<:AbstractAgentM}

w = World(myagents,myspace,p,0.)
@test size(w) ≈ 10
newa = give_birth(1,w)
addAgent!(w,newa)
@test size(w) ≈ 11
removeAgent!(w,11)
@test size(w) ≈ 10
@test isnothing(update_clock!(w,.1))

using BenchmarkTools
if false
    # TODO
    # Here we observe that deleting an agent is a bit of a pain
    # Since one has to look where this agent is located in the structure
    # This should be improved at some point
    # This could be improved by maintaining a list of agents IDs
    dicttest = Dict((i,w.agents[i]) for i in 1:p["NMax"]);
    @btime w.agents[ABMEv._get_idx(w)[8]]
    @btime dicttest[8]
    @btime ABMEv.agents(w)
end
##############  Testing world with Gillepsie
@test !isnothing(updateWorld!(w,Gillepsie()))
