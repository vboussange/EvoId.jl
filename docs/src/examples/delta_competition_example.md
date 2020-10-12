# A first model of meta population

In this script, we model the evolution of a population where agents are simply defined by their position on some landscape. We implement the simplest possible birth and death function.

## The landscape
Let's start by a linear landscape. We define a discrete segment of length `9`, with reflecting boundary conditions. In fact, reflecting boundary conditions are implemented for any finite space.

!!! warning "1D Reflections"
  As of v1, only reflections on one dimensional space are implemented. We have a two dimensional reflection method that will be released in the future.

There are two ways of implementing a linear landscape. The first one uses a `DiscreteSegment` while the second relies on LightGraphs.jl library. Both are *almost* equivalent.

### DiscreteSegment
```julia
using ABMEv
nodes = 10
mysegment = DiscreteSegment(1,nodes)
```

### grid
```julia
using ABMEv, LightGraphs
nodes = 10
g = grid([nodes,1])
mysegmentgraph = GraphSpace(g)
```
Here is how you can visualise the landscape.

```julia
using GraphPlots
gplot(g,locs_x = 1:nodes,locs_y=1:nodes)
```

## Defining competition processes
We propose that any individual have a constant birth rate, and competes with all the individuals present in the same patch. Let ``i \in \N``,``x_i \in \{1,2,\dots,9\}``.
The competition pressure experience by individual ``i`` is such that

```math
d(x_i) = \sum_j \delta(x_i-x_j)
```
where ``\delta`` is the dirac function.

In this way, we recover a logistic growth function for subpopulation within a patch.

```julia
K0 = 1000 # We will have in total 1000 individuals
b(X,t) = 1 / nodes
d(X,Y,t) = (X[1] ≈ Y[1]) / K0
```
At equilibrium, population size in each deme will reach `K0 / nodes`.

## Dispersal
We assume that anytime an offspring is born, it is given a chance to move (`\mu = 1`).
```julia
mu = [1.]
D = (1.5,)
```
## Running the world
We initialise the world with initial population of size ``K_0 / 9`` located on patch `5`. We keep track of individuals' ancestors by setting `ancestors=true`. Because we wish to use `Gillepsie` algorithm, we need `rates=true` as agents' internal birth and death rates are updated at every time step.
!!! note "Warning"
  `rates` treatment is something we might implement in the library internals.

```julia
NMax = 2000
tend = 2
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax
myagents = [Agent(myspace,(5,),ancestors=true,rates=true) for i in 1:K0/nodes]
w0 = World(myagents,myspace,p,0.)
@time sim = run!(w0,Gillepsie(),tend)
```
## Analysis
Let's verify that the population's growth is logistic. We will plot the population size over time.
To do so, one need to define `dt_saving`.

```julia
using Plots
@time sim = run!(w0,Gillepsie(),tend,dt_saving=1)
Plots.plot()
@animate for (i,t) in tspan(sim)
  Plots.
end
```

One might be tempted to plot the world for some time steps. This requires some knowledge of the library Plots.jl.
Agents values are accessed by `get_x(w::World)`.

!!! notes "Plotting"
  This will be proposed as a blackbox with PlotRecipes.jl in the near future
