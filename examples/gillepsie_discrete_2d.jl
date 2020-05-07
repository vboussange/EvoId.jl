using ABMEv,UnPack
cd(@__DIR__)
using Dates
using JLD2
using Random;Random.seed!(0)


a = 1;
D = [1e0;a];
mu = [1.;.001];
sigma_K = 1.;
sigma_a = .8;
nodes = 10
K0 = 1000
K(X) = gaussian(X[2],X[1] * a,sigma_K) / nodes
alpha(X,Y) = (X[1] ≈ Y[1]) * gaussian(X[2],Y[2],sigma_a) / K0
NMax = 5000
tend = 1.
reflected = true;
dt_saving = 15.

p = Dict{String,Any}()
@pack! p = alpha,K,D,NMax,tend,reflected,nodes, dt_saving,mu


agent0 = [Agent{MixedAgent}( Float32[rand(1:p["nodes"]), 1e-2* randn() + (p["nodes"] + 1 ) * a/2 ] ) for i in 1:K0]
world0 = vcat(agent0[:],repeat([missing],Int(p["NMax"] - K0)))
worldall,p["tspan"] = runWorld_store_G(p,world0);
# @save "worldall_discrete_a=$(a)_mu2=$(mu[2])_D2=$(D[2])_D1=$(D[1]).jld2"

using Plots
Plots.plot(worldall,p,trait=2)
