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
b(X,t) = 1.
d(X,Y,t) = gaussian(X[1],Y[1],sigma_a)/gaussian(X[1],0.,sigma_K)/K0
D = (1e-2,)
mu = [.1]
NMax = 2000
tend = 2000
# Cbar = b([0],0.)/K0 + d([0],[0],0.)
dm = d([0],[0],0.);bm = 1.
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax,dm,bm

myagents = [Agent(myspace,(-.01 + 1e-2 * randn(),),ancestors=false,rates=false) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
w1 = copy(w0)
@info "Running simulation with CFM algorithm"
@time sim = run!(w1,CFM(),tend,dt_saving=10.)

# ABMEv.clean!(sim)
using Plots
Plots.plot(sim)
