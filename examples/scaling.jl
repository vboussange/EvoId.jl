using ABMEv

## Defining parameters

sigma_K = 1.; #bandwith of resource
sigma_a = 1.2; #bandwith of competition
K0 = 1000 #carrying capacity
b(X) = gaussian(X[1],0,Float32(sigma_K)) #birth
d(X,Y) = gaussian(X[1],Y[1],Float32(sigma_a)) / Float32(K0) #competition
D = [1e-2] #mutation range
mu = [1.] #probability of mutation
NMax = 2000 #number of individual
dt_saving = 1.0 #time step saving
tend = 2.
using UnPack
p = Dict{String,Any}()
@pack! p = K,alpha,D,mu,NMax,dt_saving,tend

## Initial conditions
agent0 = [Agent(- Float32(1.5) .+ .1 .* randn(Float32,1)) for i in 1:K0]
world0 = vcat(agent0[:],repeat([missing],Int(p["NMax"] - K0)))

## launch simulation
worldall,p["tspan"] = runWorld_store_G(p,world0);

using Plots
Plots.plot(worldall,p)
