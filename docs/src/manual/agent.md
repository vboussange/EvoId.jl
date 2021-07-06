# Agent properties

## The `Agent` structure
`Agent` is the atomic structure of EvoId.jl. It has four attributes
-  the ancestors' history of traits, and the corresponding time where the traits have changed,
- a death rate and a birth rate.
```julia
mutable struct Agent{A<:Ancestors,R<:Rates,T<:Tuple,U,V} <: AbstractAgent{A,R}
    # history of traits for geotraits
    x_history::Array{T,1}
    # birth time of ancestors
    t_history::Array{U,1}
    # death rate
    d::V
    #birth rate
    b::V
end
```
!!! note "Specificities"
    The type `Agent` has two important composite types

    - `Ancestors{bool}` : when `bool = true`, the ancestors traits are stored,
    - `Rates{bool}` : when `bool = true`, the rates `d` and `b` of agents are updated at each time step. This is needed in e.g. Gillepsie Algorithm

```@autodocs
Modules = [EvoId]
Pages   = ["EvoId_Agent.jl"]
```
