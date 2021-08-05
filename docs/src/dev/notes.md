# Developping the code
I recommend to first clone your branch in the directory you like best, and then to develop, you can
```julia
using Pkg
Pkg.dev("path_to_EvoId_dir")
```
You can also do the same trick with directly the gitlab address, cf [Pkg.jl](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html)

## Future directions
- Fix `WF` algorithm
- Implement Moran process
- extend the birth function to the form `b(X,Y,t)` for coherence with death function 
- make mutation and disperal range features of agents, so that they can also evolve.
- Simplify composite type `Rates{}` parameter for `Agents`.
```julia
abstract type AbstractAgent{A<:Ancestors,R<:Rates} end
```
