# ABMEv.jl Documentation

This is a suite for simulating an Agent Based Model that captures the evolutionary dynamics of a population in a multidimensional space.

## Getting started
```julia
using ABMEv
```

## The `Agent` structure
This package is an Agent Based Model, where the atomic structure is the `Agent`. It has four attributes
-  the ancestors' history of traits, and the corresponding time where the traits have changed,
- a death rate and a birth rate.
```julia
mutable struct Agent{T,U}
    # history of traits for geotraits
    x_history::Array{U}
    # birth time of ancestors
    t_history::Array{U,1}
    # death rate
    d::Float64
    #birth rate
    b::Float64
end
```
For more understanding of the composite types `T,U`, check the wiki: Gillepsie / Agent Type.
## Parameters of the simulation
Parameters are stored in the parameter dictionary `p`
### General parameters
- ```"reflected"=>true```: if ```true``` then reflection occurs on the first trait -which should stand for geographic position. Depending on the agent type, reflections occurs in the domain $` [-1,1] `$ or between nodes 1 and `p["nodes"]`
- ```"alpha" => α```: is the competition function
- ```"K" => K```: is the birth rate
- ```"tend" => 1.5```: is the time to end simulation

:warning: Check how to define functions α and K in the algorithm section.
### Mutation
If anisotropy in mutation, the following parameters should be declared as arrays where each entry corresponds to a dimension.
- ```mu``` The probability of mutation.
- ```D``` If mutation happens on the agent, the increment follows $`\mathcal{N}_{ 0, D}`$
### Birth
#### Growth
- ```K``` is the birth coefficient ( $`b(x) = K(x)`$ )
### Death
#### Competition
- Competition between agent with trait ```x``` and ```y``` is defined as
```α(x,y)```

- Death coefficient is defined as $`d(x^{(i)}) = \sum_j^{N(t)} \alpha(x^{(i)},x^{(j)})`$

### Fitness
Fitness is defined internally as ```b - d```.
> TODO ```b``` here is confounded with ```K```.


## Launching simulation
Two type of simulation algorithm can be used
```@contents
Pages = ["manual/gillepsie.md",
        "manual/wright_fisher.md"]
Depth = 5
```
### Gillepsie algorithm

```julia
using ABMEv
## Defining parameters
sigma_K = 1.; #bandwith of resource
sigma_a = 1.2; #bandwith of competition
K0 = 1000 #carrying capacity
K(X) = gaussian(X[1],0,Float32(sigma_K)) #birth
alpha(X,Y) = gaussian(X[1],Y[1],Float32(sigma_a)) / Float32(K0) #competition
D = [1e-2] #mutation range
mu = [1.] #probability of mutation
NMax = 2000 #number of individual
dt_saving = 1.0 #time step saving
tend = 1000.
using UnPack
p = Dict{String,Any}()
@pack! p = K,alpha,D,mu,NMax,dt_saving,tend

## Initial conditions
agent0 = [Agent(.1 .* randn(Float32,1)) for i in 1:K0]
world0 = vcat(agent0[:],repeat([missing],Int(p["NMax"] - K0)))

## launch simulation
worldall,p["tspan"] = runWorld_store_G(p,world0);

using Plots
Plots.plot(worldall,p)
```


### Wright Fisher algorithm
```julia
sigma_K = .9;
sigma_a = 1.251;
K0 = 1000;
# K(X) = gaussian(X[1],0.,sigma_K)
K(X) = 1 - 0.125 * sum(X.^2)
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
# α(X,Y) = 0.
p = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 10.)
na_init = K0
agents = [Agent( [1e-2]  .* randn(1) .- .5) for i in 1:K0]
@time worldall_test,p["tspan"] = runWorld_store_WF(p,agents,mode="std");
```




## References
- [Champagnat and Ferriere founding article](https://linkinghub.elsevier.com/retrieve/pii/S0040580905001632)
- [Champagnat and Ferriere second article - 2008](https://www.tandfonline.com/doi/full/10.1080/15326340802437710)

## Similar packages:
[Agents.jl](https://juliadynamics.github.io/Agents.jl/) It would be worth to have a look! It has been designed by Ali R. Vahdati, from UZH.
