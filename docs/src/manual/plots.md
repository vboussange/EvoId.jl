# Plotting

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
- ```"gs"``` returns a scatter plot `(xaxis = geotrait, yaxis = trait value)`
- ```"3dgeo"``` plots a 3d diagram with x axis as geotrait and y axis as the second component. :warning: this is probably weekend
- ```"3d"``` plots a 3d diagram with first and second component as x and y axis
- ```"var"``` plots the variance of the  component specified by ```trait=2``` :question: with respect to time?
- ```"vargeo"``` plots the variance of the geotrait

```@autodocs
Modules = [ABMEv]
Pages   = ["ABMEv_plotting.jl"]
```
