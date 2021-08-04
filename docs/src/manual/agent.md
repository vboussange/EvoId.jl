# Agent properties

!!! note "Specificities"
    The type `Agent` has two important composite types

    - `Ancestors{bool}` : when `bool = true`, the ancestors traits are stored,
    - `Rates{bool}` : when `bool = true`, the rates `d` and `b` of agents are updated at each time step. This is needed in e.g. Gillepsie Algorithm

```@autodocs
Modules = [EvoId]
Pages   = ["Agent.jl"]
```
