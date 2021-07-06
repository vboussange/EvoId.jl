# Wright  Fisher algorithm
## Foundations

The Wright Fisher process is an individual based model where the number of agents is constant through time. It is helpful to visualize it through marbles in jar:

![alt text](https://upload.wikimedia.org/wikipedia/commons/0/0b/Random_sampling_genetic_drift.svg)

At each time step, ``N`` agents are picked up from previous generation to reproduce. Their number of offspring is proportional to their fitness, calculated as usual with **birth and death rates**.

It takes thus **only one time step to go trough one generation**. Thus it is more suit- able for numerical simulations. In practice, the Moran and Wright–Fisher models give qualitatively similar results, but genetic drift runs twice as fast in the Moran model.


From this perspective we can easily get that branches are less stable than in the Gillepsie scenario, for as as time goes to infinity the probability of going extinct is intuitively bigger than 0.


## Initial conditions
One need to construct the world as an array of agents, which will be the ancestors of the following
```julia
    agents = [Agent( [2e-1]  .* randn(1)) for i in 1:K0]
```
The function

is then called `p["tend"] -1` times.
## Scenarios
You have several options available concerning the resource implemented and competition:
- ``` mode="std"``` is the standard mode
- ``` mode="grad2D"``` corresponds to a an ecological gradient
>We are not sure whether this corresponds to the following two images
- ``` mode="mountain"``` corresponds to a scenario where a mountain arises (with an ecological gradient)
- ``` mode="split"``` corresponds to a scenario where the resource is splitted in two
- ``` mode="graph"``` this guy is probably not working

## Parallelism
You can run your script in parallel, which makes sense for large populations. To do so:
```julia
using Distributed;addprocs()
@everywhere using EVOID
```
Parallelism only works with Wright Fisher model.


```@autodocs
Modules = [EVOID]
Pages   = ["EVOID_WF.jl"]
```
