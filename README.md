# ABMEv.jl Documentation
This is a suite for simulating an Agent Based Model that captures the evolutionary dynamics of a population in a multidimensional space.

## Installation
```julia
using Pkg;
Pkg.add("https://gitlab.ethz.ch/bvictor/abmev.git")
```
This will download latest version from git repo and download all dependencies.
To check out from an other branch than master, one has to do the trick
```julia
using Pkg;
Pkg.add("ABMEv#no_C_matrix")
```
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
As of now, no mode is implemented. For further examples, check the folder `examples` in source code.
#### Specific parameters

- ```dt_saving = 10.```
will allow to save the world every 10. time steps. If not specified, the algorithm will return first and last time step world.
- ```NMax``` Maximum number of individuals that can be attained. If attained, then the programm stops.

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
You have several options available concerning the resource implemented and competition:
- ``` mode="std"``` is the standard mode
- ``` mode="grad2D"``` corresponds to a an ecological gradient
>We are not sure whether this corresponds to the following two images
- ``` mode="mountain"``` corresponds to a scenario where a mountain arises (with an ecological gradient)
- ``` mode="split"``` corresponds to a scenario where the resource is splitted in two
- ``` mode="graph"``` this guy is probably not working


### Parallelism
You can run your script in parallel, which makes sense for large populations. To do so:
```julia
using Distributed;addprocs()
@everywhere using ABMEv
```
Parallelism only works with Wright Fisher model.
## Properties of agents
### Main accessors
You can access properties of the agent using the following functions

- `get_geo(a::Agent{U,T},t::Number) where {U,T}`: returns the geotrait
- `get_x(a::Agent,t::Number,i::Integer)` returns trait `i` of the agent. Geotrait corresponds to dimension `i=0`. If you do not want to access geotrait, you can also use `get_x(a::Agent,i::Integer)`.

### World accessors
You can access properties of the world of agents using the following functions

- `get_x(world::Array{T},t::Number,trait::Integer) where {T <: Agent}`

Returns trait of every agents of world in the form of an array which dimensions corresponds to the input. If trait = 0 , we return the geotrait. One can also use `get_x(world::Array{T},trait::Integer) where {T <: Agent}` if no need for geotrait
- `get_xarray(world::Array{T,1},t::Number,geotrait::Bool=false) where {T <: Agent}`

Returns every traits of every agents of `world` in the form **of a one dimensional array** (in contrast to `get_x`). If `geotrait=true` the geotrait is also added to the set of trait, in the last line. If you do not want to specify `t` (only useful for geotrait), it is also possible to use `get_xarray(world::Array{T,1}) where {T <: Agent}`.


Returns every traits of every agents of world in the form of an array
If geotrait = true, then a last trait dimension is added, corresponding to geotrait.

- `get_xhist(world::Vector{Agent},geotrait = false)` : :warning: This method is broken.

Returns the trait history of every agents of world in the form of an 3 dimensional array,
with

    - first dimension as the agent index
    - second as time index
    - third as trait index
If geotrait = true, then a last trait dimension is added, corresponding to geotrait.
Note that because number of ancestors are different between agents, we return an array which size corresponds to the minimum of agents ancestors,
and return the last generations, dropping the youngest ones

- `world2df(world::Array{T,1},t::Number,geotrait = false) where {T <: Agent}`

Converts the array of agent world to a dataframe, where each column corresponds to a trait of the agent, and an extra column captures fitness.
Each row corresponds to an agent.
Again, if one do not want to feed `t` (no need for geotrait), it is possible to use `world2df(world::Array{T,1}) where {T <: Agent}`.

### Other accessors
> TODO: describe the following accessors
```julia
get_x(a::Agent,i::Number) = a.x_history[Int(i),end]
get_x(a::Agent) = a.x_history[:,end]
get_xhist(a::Agent,i::Number) = a.x_history[Int(i),:]
get_xhist(a::Agent) = a.x_history
get_geo(a::Agent) = sum(get_xhist(a,1))
get_d(a::Agent) = a.d
get_b(a::Agent) = a.b
get_fitness(a::Agent) = a.b - a.d
get_dim(a::Agent) = size(a.x_history,1)
get_nancestors(a::Agent) = size(a.x_history,2)
```
> Todo : implement geotrait in `get_xhist`

## Plotting
ABMEv comes with Plot recipes:
`plot(world::Array{U},p;what=["x","H"],trait = 1) where U <: Union{Missing,Agent{T}} where T`.

An example to use it: 
```julia
using Plots;pyplot()
Plots.plot(worldall,p_default,what=["x"],trait=2)
```
You can specify what you want to plot in the array ```what```:
- ```"x"``` returns a scatter plot `(xaxis = time, yaxis = trait)` . Geotrait corresponds to `trait=0`
- `"xs"` only works for agent type `MixedAgent`, because it needs a discrete geographical space. It returns a scatter plot `(xaxis = geographical component, yaxis = trait value)`. Very similar to a `histogram2d` plot, with nicer look.

- ```"3dgeo"``` plots a 3d diagram with x axis as geotrait and y axis as the second component. :warning: this is probably weekend
- ```"3d"``` plots a 3d diagram with first and second component as x and y axis
- ```"var"``` plots the variance of the  component specified by ```trait=2``` :question: with respect to time?
- ```"vargeo"``` plots the variance of the geotrait 

## Developping the code
I recommend to first clone your branch in the directory you like best, and then to 
To develop, you ca
```julia
using Pkg
Pkg.dev("path_to_ABMEv_dir")
```
You can also do the same trick with directly the gitlab address, cf [Pkg.jl](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html)

## References
- [Champagnat and Ferriere founding article](https://linkinghub.elsevier.com/retrieve/pii/S0040580905001632)
- [Champagnat and Ferriere second article - 2008](https://www.tandfonline.com/doi/full/10.1080/15326340802437710)
