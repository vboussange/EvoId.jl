p = Dict(1:1:1000 .=> 1:1:1000)
delete!(p,999)
a = Any[rand() for i in 1:1000]
a[2:2:1000] .= missing
get_idx(w) = collect(eachindex(skipmissing(w)))
using BenchmarkTools
@btime a[get_idx(a)[3]]
@btime p[3]


@btime [true for i in 1:1000]
[true for i in 1:1000]

@btime sort!(values(p))
function findfreeidx(p)
    i = 1
    while haskey(p,i)
        i+=1
    end
    i
end

function get_xp(a,j)
    for (i,a) in enumerate(skipmissing(a))
        if i==j
            return a
            break
        end
    end
end
@btime get_xp(a,998)

@btime findfreeidx(p)
@btime collect(keys(p))[2]

@btime(count(ismissing,))

############################################################
###########Reflection on vectors of missing or not##########
############################################################
cd(@__DIR__)
using Random
Random.seed!(0)
using LightGraphs
using Test
using Revise,ABMEv
using UnPack,JLD2

myspace = (RealSpace{1,Float64}(),)
sigma_K = .9;
sigma_a = .7;
K0 = 1000;
b(X) = gaussian(X[1],0.,sigma_K)
d(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
D = (1e-2,)
mu = [.1]
NMax = 10000000
tend = 1.5
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax

using BenchmarkTools
a = Agent(myspace,(0,),ancestors=true,rates=true)
myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
@btime push!(myagents,a)

###### Testing world adding agents
myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
size(w0)
@btime addAgent!(w0,a)
size(w0)
@btime removeAgent!(w0,1)

###### Testing sim adding agents
using BenchmarkTools

myspace = (GraphSpace(SimpleGraph(10,10)),RealSpace{1,Float64}())
myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
w1 = World(myagents,myspace,p,0.)
@time sim = run!(w1,Gillepsie(),tend,dt_saving = .1)
get_size(sim)
# @btime add_entry!(sim,w1)
