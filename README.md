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
Two type of simulation algorithm can be used

### Gillepsie algorithm

```julia
a = 0;
sigma_K = .9;
sigma_a = .7;
K0 = 1000;
K(X) = gaussian(X[1],0.,sigma_K)
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
p_default = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 1.5,
        "NMax" => Int(10000))
na_init = K0
world0 = new_world_G(na_init,p_default,spread = [.01 .01], offset = [-.5 -.5])
worldall,p_default["tspan"] = runWorld_store_G(p_default,world0)
```
As of now, no mode is implemented.
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

## Developping the code
I recommend to first clone your branch in the directory you like best, and then to 
To develop, you ca
```julia
using Pkg
Pkg.dev("path_to_ABMEv_dir")
```
You can also do the same trick with directly the gitlab address, cf [https://docs.julialang.org/en/v1/stdlib/Pkg/index.html](Pkg.jl)

## How does it implement mechanisms ?

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

## References
- [Champagnat and Ferriere founding article](https://linkinghub.elsevier.com/retrieve/pii/S0040580905001632)
- [Champagnat and Ferriere second article - 2008](https://www.tandfonline.com/doi/full/10.1080/15326340802437710)
