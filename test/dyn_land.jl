cd(@__DIR__)
using Random
Random.seed!(0)
using LightGraphs
using Test
using Revise,ABMEv
using UnPack,JLD2

g1 = LightGraphs.grid([9,1])
g2 = SimpleGraph(9)
# This is the function to implement DynGraphSpace.
# Note that it returns indices
g = [g1,g2]
function periodic_update_graph(T,t)
    if sin(t/T*2*π) > 0
        1
    else
        2
    end
end
update_g(t) = periodic_update_graph(10,t)
dyng = DynGraphSpace(g,update_g)
## testing atomic methods
@test  last(randomwalk(ABMEv.get_graph(dyng,16.),1,10)) ≈ 1
@test  !(last(randomwalk(ABMEv.get_graph(dyng,1.),1,10)) ≈ 1)

## simulations
myspace = (dyng,)
sigma_K = .9;
sigma_a = .7;
K0 = 1000;
b(X,t) = gaussian(X[1],0.,sigma_K)
d(X,Y,t) = gaussian(X[1],Y[1],sigma_a)/K0
D = (1e-2,)
mu = [.1]
NMax = 10000
tend = 1.5
p = Dict{String,Any}();@pack! p = D,mu,NMax
myagents = [Agent(myspace,(1,),ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
w1 = copy(w0)
@info "Running simulation with Gillepsie algorithm"
sim = run!(w1,Gillepsie(),tend,b,d)

@test !(isnothing(sim))
