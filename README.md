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
To develop, you can add the dependencies in the project by doing
```julia
using Pkg
Pkg.activate("path_to_ABMEv")
```

## Birth and Death mechanisms
> We are always balanced between taking the integral of the competition and resource kernel as constant, or taking its maximum peak as constant.

:poop:
## Geotrait
The geotrait is calculated *a posteriori*, and is not taken into account during the simulation.
> It used to be but for the sake of simplicity we now forget about it.

### Mutation
If anisotropy in mutation, the following parameters should be declared as arrays where each entry corresponds to a dimension.
- ```mu``` The probability of mutation.
- ```D``` If mutation happens on the agent, the increment follows $\mathcal{N}_{ 0, D}$
### Birth
#### Growth
- Resource kernel for agent with trait $x$ is defined as 
```math
K_{\mu,\sigma}(x) = K_0 \mathcal{N}_{\mu,\sigma}(x)
```
with $\mu$ and $\sigma$ potentially vectors.
> We just modified this in ABMEv_Agent.jl so you should check if it works.
- Dirth coefficient is defined as $`b(x) = K(x)`$
### Death
#### Competition
- Competition between agent with trait ```x``` and ```y``` is defined as
```math 
\alpha(x,y) = \exp(-\sum_i^{N(t)} \frac{1}{\sigma_{\alpha_i}^{n_\alpha}}\sum_j^T (x_{i,j} - y_{i,j})^{n_\alpha})
```
- Death coefficient is defined as $`d(x^{(i)}) = \sum_j^{N(t)} \alpha(x^{(i)},x^{(j)})`$
> We are not sure if the sum includes $`x^{(i)}`$ or not.
### Fitness
Fitness is defined as ```b - d```.

## Parameter description
- ```K0``` Carrying capacity
- ```a``` only used for mode ```grad2D``` where growth rate is set such as $`\mu = a x_1`$
> We are not sure if this is OK or not? Check it 
[Grad2D kernel explained](https://gitlab.ethz.ch/bvictor/abmev/-/wikis/Grad2D)
### Gillepsie algorithm

```julia
p_default = Dict("K0" => 1000.,
        "D" => [1e-2 - 1e-3],
        "mu" => [1. .1],
        "sigma_a" => [5e-1; 7e-1],
        "sigma_K" => [9e-1; 9e-1],
        "n_alpha" => 2.,
        "n_K" => 2.,
        "tend" => 150.,
        "NMax" => Int(2e3))
na_init = 500
world0 = new_world_G(na_init,p_default,spread = [.01 .01], offset = [-.5 -.5])
worldall,p_default["tspan"] = runWorld_store_G(p_default,world0)
```
As of now, no mode is implemented.

### Wright Fisher algorithm
#### Specific parameters
- ```NMax``` Maximum number of individuals that can be attained. If attained, then the programm stops.
```julia
p_default = Dict("K0" => 2000.,
        "D" => [2e-2], # this is the dispersion coefficient, that should be taken as constant
        "mu" => [1.],
        "sigma_a" => [2e0],
        "sigma_K" => [1e0],
        "a" => .95,
        "n_alpha" => 2.,
        "n_K" => 2.,
        "tend" => 10.,
        "NMax" => Int(1e4))
na_init = p_default["K0"]
agent0 = [Agent([-.5 + .01  * randn()]) for i in 1:na_init]
world0 = vcat(agent0[:],repeat([missing],Int(p_default["NMax"] - na_init)))
(worldall , p_default["tspan"]) = runWorld_store_WF(p_default,agent0,mode="std")
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
@everywhere push!(LOAD_PATH,homedir()*"path_to_ABMEv")
@everywhere using ABMEv
```

## Properties of agents
You can access properties of the agent using the following functions
- `get_xarray(world::Array{Agent{T}},trait::Int) where T`

Returns trait of every agents of world in the form of an array

> TODO: describe the following accessors
```julia
get_x(a::Agent,i::Number) = a.x_history[Int(i):Int(i),end]
get_x(a::Agent) = a.x_history[:,end]
get_xhist(a::Agent,i::Number) = a.x_history[Int(i):Int(i),:]
get_xhist(a::Agent) = a.x_history
get_geo(a::Agent) = sum(get_xhist(a,1))
get_d(a::Agent) = a.d
get_b(a::Agent) = a.b
get_fitness(a::Agent) = a.b - a.d
```

## Plotting
ABMEv comes with Plot recipes:
`function plot(world::Array{U},p;what=["x","H"],trait = 1) where U <: Union{Missing,Agent{T}} where T`.

An example to use it: 
```julia
using Plots;pyplot()
Plots.plot(worldall,p_default,what=["x"],trait=2)
```
You can specify what you want to plot in the array ```what```:
- ```x``` plots the component specified by ```trait=2```
- ```geo``` plots geotrait, computed from first component
- ```3dgeo``` plots a 3d diagram with x axis as geotrait and y axis as the second component
- ```3d``` plots a 3d diagram with first and second component as x and y axis
- ```var``` plots the variance of the  component specified by ```trait=2```
- ```vargeo``` plots the variance of the geotrait
## References
- [Champagnat and Ferriere founding article](https://linkinghub.elsevier.com/retrieve/pii/S0040580905001632)
- [Champagnat and Ferriere second article - 2008](https://www.tandfonline.com/doi/full/10.1080/15326340802437710)
