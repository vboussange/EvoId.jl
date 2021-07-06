# Developping the code
I recommend to first clone your branch in the directory you like best, and then to
To develop, you ca
```julia
using Pkg
Pkg.dev("path_to_EVOID_dir")
```
You can also do the same trick with directly the gitlab address, cf [Pkg.jl](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html)

## Future directions
- Try to improve parallelism with the help of Threads.@Threads and @inbounds (cf [tutorial](https://juliagpu.gitlab.io/CUDA.jl/tutorials/introduction/#Introduction-1) )
- Make use of CUDA.jl to accelerate the simulations wih the use of GPU.

## Todo
```julia
abstract type AbstractAgent{A<:Ancestors,R<:Rates} end
```
It seems that we do not need to have the `Rates{}` parameter for `Agent` type.


## Notes
### Numerics
:warning:
[Don’t pass expressions, or strings, pass functions](https://www.youtube.com/watch?v=mSgXWpvQEHE&t=573s)

### Birth and Death mechanisms
<!-- > We are always balanced between taking the integral of the competition and resource kernel as constant, or taking its maximum peak as constant. -->

<!-- :poop: -->
### Geotrait
The geotrait is calculated *a posteriori*, and is not taken into account during the simulation.
> It used to be but for the sake of simplicity we now forget about it.
