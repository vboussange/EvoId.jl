# ABMEv.jl Documentation

This is a suite for simulating the evolutionary dynamics of a population in a multidimensional space. The population is modelled at the level of the individual - this suite hence falls in the realm of *Agent Based Modelling*.
The purpose of this package is to provide a numerical laboratory for evolutionary dynamics, supplying
- a flexible atomic structure of agents and underlying evolutionary space
- algorithms and update rules for the simulations
- analysis tools to investigate the simulations.

## Features
Agents consist of sets of traits in some combination of vector spaces. A vector space can represent for example **a geographical landscape, a trait space, or genetic structure**. Spaces can be of any dimensions, discrete or continuous, bounded or unbounded. They can equally consist of graphs.
Vector spaces are used to define birth and death processes, as well as mutation processes.

### Specificities
- [ABMEv.jl allows to keep track of agents' trait lineages](@ref lineages)
- [ABMEv.jl enables to run evolutionary dynamics on graphs!](@ref genetic_structure)

## Getting started
```@repl
using ABMEv
```

## Tutorial
We strongly advise to have a look at the tutorial section.
```@contents
Pages = [
    "examples/delta_competition_example.md",
    "examples/changing_environment.md",
    "examples/sympatric_speciation.md",
    "examples/gradient.md",
    "examples/genetic_structure.md",
    ""
    ]
Depth = 2
```
## How it works
General workflow to launch any simulation is the following

- [Define the combination of vector spaces you are interested in.](manual/space.md)
- Define birth and death function, that depend on agents position in the space
- Define mutation function
- [Define initial population state and time](manual/world)
- [Run the simulation according to some updating algorithm](manual/run_world.md)
- [Obtain a summary of the population state](manual/callbacks.md)

### Available algorithms
As of now, three types of simulation algorithm can be used:
- [Gillepsie](manual/gillepsie.md)
- [Wright-Fisher](manual/wright_fisher.md)
- [CFM](CFM.md)

## References
- [Champagnat and Ferriere founding article](https://linkinghub.elsevier.com/retrieve/pii/S0040580905001632)
- [Champagnat and Ferriere second article - 2008](https://www.tandfonline.com/doi/full/10.1080/15326340802437710)

## Similar packages:
[Agents.jl](https://juliadynamics.github.io/Agents.jl/) This package is oriented towards general ABM modelling, and miss many of the specificities of ABMEv.jl , such as vector spaces operations, update algorithm, and analysis tools.
