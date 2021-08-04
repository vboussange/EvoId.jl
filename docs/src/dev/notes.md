# Developping the code
I recommend to first clone your branch in the directory you like best, and then to
To develop, you can
```julia
using Pkg
Pkg.dev("path_to_EvoId_dir")
```
You can also do the same trick with directly the gitlab address, cf [Pkg.jl](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html)

## Future directions
- Try to improve parallelism with the help of Threads.@Threads and @inbounds (cf [tutorial](https://juliagpu.gitlab.io/CUDA.jl/tutorials/introduction/#Introduction-1) )
- Make use of CUDA.jl to accelerate the simulations wih the use of GPU.

## Todo
- We do not need to have the `Rates{}` parameter for `Agent` type.
```julia
abstract type AbstractAgent{A<:Ancestors,R<:Rates} end
```
- extend the birth function to the form `b(X,Y,t)` for coherence with death function 
- make mutation and disperal range features of agents, so that they can also evolve.
