# ABMEv.jl
<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://vboussange.github.io/ABMEv.jl/stable) -->
<!-- For now we only direct to dev documentation. In the future, one will need to deploy a ssh key to and use TagBot. -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://vboussange.github.io/ABMEv.jl/dev)
[![pipeline status](https://gitlab.ethz.ch/bvictor/abmev/badges/master/pipeline.svg)](https://gitlab.ethz.ch/bvictor/abmev/-/commits/master)

<div align="center"> <img
src="https://vboussange.github.io/images/research/conceptual_onlyadapt.png"
alt="" width="400"></img> </div>

ABMEv.jl is a package aimed at simulating the eco-evolutionary dynamics of a population in a multidimensional space, at the individual level. 

Individuals are characterised by a set of traits in some combination of evolutionary spaces. An evolutionary space can represent for example a geographical landscape, a trait space, or genetic structure. Individuals experience three elementary events : **birth, death, mutation (or migration)**.

The user can provide **any birth and death functions**, which should depend on the system state and the individuals' trait. Together with the mutation process, this defines the dynamics of the system.

ABMEv.jl provides a **numerical laboratory** for eco-evolutionary dynamics, supplying
- flexible types for **individuals**, which can
    - evolve over any combination of space,
    - store ancestors trait,
- flexible types for **evolutionary spaces**, that can consist of multidimensional **discrete or continuous domains**, as well as **graphs**,
- the possibility to use **callback functions** to save the state of the system at any time step
- several **algorithms** for the simulations (Gillepsie, CFM),
- **utility functions** to analyse simulation results.

## Installation
Open Julia in your favorite REPL and type the following 

```julia
using Pkg;
Pkg.add("https://gitlab.ethz.ch/bvictor/abmev.git")
```

This will download latest version from git repo and download all dependencies.

## Getting started
Check out the great documentation if you want to use the advance features of ABMEv.jl. Otherwise, you can content yourself with the simple tutorial prodived below

## A first tutorial
We provide here a tutorial a simple scenario, where the population is structured over the vertices of a graph, where each vertex selects for an optimal trait value $`\theta_i \in \{-1,1\}`$.

Let's import ABMEv.jl, and LightGraphs.jl
```julia
using ABMEv
using LightGraphs.jl
```

Now define the trait space, and the graph (i.e. the abstraction of the landscape). We choose a star graph with 7 vertices.

```julia
nodes = 7
g = star_graph(nodes)
landscape = GraphSpace(g)
θ = [rand([-1,1]) for i in 1:nodes]
traitspace = RealSpace(1)
evolspace = (landscape,traitspace)
```

Now let's define the birth and death functions. Individuals are characterised by their position on the graph and the adaptive trait.

```julia
K = 1000
b(X,t) = 1 - 0.5 * (θ[X[1]] - X[1])^2
d(X,Y,t) = (X[1] ≈ Y[1]) / K
```

Define the mutation processes. Assuming a random walk of length 1, you should specify `nothing`.
```julia
D = [nothing,5e-1]
mu = [1.,1.]
```

Simulation time, and callback function

```
tend = 200
t_saving_cb = collect(range(0.,tend,length=300))
cb() = Dict("N" => size(w))
```

Initial conditions

```julia
myagents = [Agent(evolspace,[rand(1:nodes),randn() * D[2]]) for i in 1:K]
```

And off we go

```julia
w0 = World(myagents,evolspace,p,0.)
sim = run!(w0,Gillepsie(),tend,b,d,dt_saving=2.)
```

Let's plot some cool stuff. 
Population size
<div align="center"> <img
src="https://github.com/vboussange/vboussange.github.io/blob/master/images/software/PDE_vs_ABM.png?raw=true"
alt="" width="400"></img> </div>

Population trait over time
