# ABMEv.jl Documentation
This is a suite for simulating an Agent Based Model that captures the evolutionary dynamics of a population in a multidimensional space.

## Getting started
```julia
using ABMEv
```
To develop, you can add the dependencies in the project by doing
```julia
using Pkg
Pkg.activate("path_to_ABMEv")
```
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

### Plotting
ABMEv comes with Plot recipes
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
