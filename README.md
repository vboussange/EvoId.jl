# EvoId.jl
<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://vboussange.github.io/EvoId.jl/stable) -->
<!-- For now we only direct to dev documentation. In the future, one will need to deploy a ssh key to and use TagBot. -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://vboussange.github.io/EvoId.jl/dev)
<!-- [![Build Status](https://github.com/vboussange/EvoId.jl/workflows/CI/badge.svg)](https://github.com/vboussange/EvoId.jl/actions) -->

<div align="center"><img src="docs/src/assets/gif_evoid.gif" width="400"></img> </div>

EvoId.jl (for **Evo**lutionary **I**n**d**ividual-based models) is a package aimed at simulating the eco-evolutionary dynamics of a population in a multidimensional space, at the individual level. The dynamics is specified under the framework of [stochastic models for structured populations](https://arxiv.org/abs/1506.04165).

Individuals are characterised by **a set of traits** in some **combination of evolutionary spaces**. An evolutionary space can represent for example a geographical landscape, a trait space, or genetic structure. Spaces can be of any dimensions, discrete or continuous, bounded or unbounded. They can equally consist of graphs. Individuals give birth at a rate given by the birth function `b`, and die at a rate given by the death function `d`. When an individual give birth, its offspring can move on the underlying evolutionary spaces. The movement can capture whether migration or mutation processes, and is characterised by a probability `m` and movement range `D`.

The user can provide **any birth and death functions**, which should depend on the system state and the individuals' trait. Together with the **movement rate and movement range**, this defines the dynamics of the system.

EvoId.jl provides a **numerical laboratory** for eco-evolutionary dynamics, supplying

- flexible types for **individuals**, which can
    - evolve over any combination of space,
    - [store ancestors trait](https://vboussange.github.io/EvoId.jl/dev/examples/gradient.html#lineages),
- flexible types for **evolutionary spaces**, that can consist of multidimensional **discrete or continuous domains**, as well as **graphs**,
- the possibility to use **callback functions** to save the state of the system at any time step
- several **algorithms** for the simulations ([Gillespie](https://en.wikipedia.org/wiki/Gillespie_algorithm),[Wright-Fisher](https://en.wikipedia.org/wiki/Moran_process), etc...),
- **utility functions** to analyse simulation results.

## Installation
Open Julia in your favorite REPL and type the following

```julia
using Pkg;
Pkg.add("https://github.com/vboussange/EvoId.jl")
```

This will download latest version from git repo and download all dependencies.

## Getting started
Check out the tutorial prodived below. You can also look at the `example` folder, or dive into the [documentation](https://vboussange.github.io/EvoId.jl/dev) if you want to use the advanced features of EvoId.jl. 

## Related papers
- [Topology and habitat assortativity drive neutral and adaptive diversification in spatial graphs](https://www.biorxiv.org/content/10.1101/2021.07.06.451404v2), Boussange et al. 2021.

## Similar packages
[Agents.jl](https://juliadynamics.github.io/Agents.jl/) is a library oriented towards general ABM modelling, and thus is not as easy to deploy as EvoId.jl for simulating stochastic models of structured populations.

## Contributing 
Please feel free to contact me! :)


-----
## Tutorial
We provide here a tutorial that sums up the 5 steps necessary to launch a simulation. For the sake of the tutorial, we propose to model a population that is structured over the vertices of a graph and characterised by a trait under selection.

### Copy paste code
```julia
using EvoId
using LightGraphs
using Plots

nodes = 7
g = star_graph(nodes)
landscape = GraphSpace(g) # spatial space
θ = [rand([-0.5, 0.5]) for i in 1:nodes] # optimal trait values
# X[1][] is the geographical position
# X[2][] corresponds to the adaptive traits
traitspace = RealSpace(1) # phenotypic space
evolspace = (landscape, traitspace) # space over which individuals are structured

const K = 100. # carrying capacity
b(X,t) = 1. - 0.5 * (θ[Int(X[1][])] - X[2][])^2 # birth function
d(X,Y,t) = (X[1] ≈ Y[1]) ? 1. / K : 0. # death function
alg = Gillepsie() # update rule
NMax = 2000 # maximum number of individuals
D = [nothing, 5e-2] # dispersion coefficient
mu = [1e-1, 1e-1] # mutation / migration rate


myagents = [Agent(evolspace, [rand(1:nodes,1), randn(1) * D[2]]) for i in 1:K] # array containing the founder individuals
# random position on the graph
# random position on the trait space centered around 0
w0 = World(myagents, evolspace, D, mu, NMax, 0.) # the initial world, defined at time 0.

tend = 500. # time horizon
t_saving_cb = collect(range(0., tend, length=200)) # time step where callback function is called
cb(w) = Dict("N" => length(w)) # callback function

println("Running simulation with callback function")
@time sim = run!(w0, alg, tend, b, d, cb=cb, t_saving_cb = t_saving_cb)
Plots.plot(sim.tspan, sim["N"])

println("Running simulation with `dt_saving`")
myagents = [Agent(evolspace, [rand(1:nodes,1),randn(1) * D[2]]) for i in 1:K]
w0 = World(myagents, evolspace, D, mu, NMax, 0.)
@time sim = run!(w0, Gillepsie(), tend, b, d, dt_saving = 2.0)
# Plots.plot(sim, trait=2)
```

### Details
#### 0. Import the relevant libraries
Let's import EvoId.jl, and LightGraphs.jl
```julia
using EvoId
```

#### 1. Define the evolutionary spaces
We define the geographical space as star graph with 7 vertices (i.e. the abstraction of the landscape), and a continuous trait space.

```julia
nodes = 7
g = star_graph(nodes)
landscape = GraphSpace(g) # spatial space
θ = [rand([-0.5, 0.5]) for i in 1:nodes]
traitspace = RealSpace(1) # phenotypic space
evolspace = (landscape, traitspace) # space over which individuals are structured
```

#### 2. Define birth and death function
Birth and death functions depend on individuals position in the combination of spaces defined above, i.e. position on the graph and the adaptive trait.
We decide that each vertex selects for an optimal trait value $`\theta_i \in \{-1,1\}`$.

```julia
const K = 100. # carrying capacity
b(X,t) = 1. - 0.5 * (θ[Int(X[1][])] - X[2][])^2 # birth function
d(X,Y,t) = (X[1] ≈ Y[1]) ? 1. / K : 0. # death function
alg = Gillepsie() # update rule
NMax = 2000 # maximum number of individuals
D = [nothing, 5e-2] # dispersion coefficient
mu = [1e-1, 1e-1] # mutation / migration rate
```
> :warning: birth and death functions should have the same number of arguments as above.

#### 3. Define how individuals move over the evolutionary space
Individual movements correspond to migration and mutation processes. On continuous spaces, one should specify a migration range and a migration rate. On discrete spaces, only a migration rate is needed (one assumes that indivuals can migrate only to neighbour patches).

```julia
NMax = 2000 # maximum number of individuals
D = [nothing, 5e-2] # dispersion coefficient
mu = [1e-1, 1e-1] # mutation / migration rate
```

#### 4. Define the initial population state

```julia
myagents = [Agent(evolspace, [rand(1:nodes,1), randn(1) * D[2]]) for i in 1:K] # array containing the founder individuals
# random position on the graph
# random position on the trait space centered around 0
w0 = World(myagents, evolspace, D, mu, NMax, 0.) # the initial world, defined at time 0.
```

#### 5. Run
Simulation time, and callback function

```julia
tend = 500. # time horizon
t_saving_cb = collect(range(0., tend, length=200)) # time step where callback function is called
cb(w) = Dict("N" => length(w)) # callback function
```


And off we go

```julia
@time sim = run!(w0, alg, tend, b, d, cb=cb, t_saving_cb = t_saving_cb)
```
#### Plotting
```julia
using Plots
plot(sim.tspan, sim["N"])
```

![](docs/src/assets/tutorials/delta_comp_wsize.png)

With a few more tricks, one can also plot the population trait density over time, for example the local trait density for individuals living on vertex 1.

![](docs/src/assets/ABM_local_trait_dens_adapt.png)

Check out the folder `examples` in the git repo to see this tutorial in a julia file, as well as plenty others!
