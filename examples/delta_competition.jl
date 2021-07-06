
# Space as a discrete segement
using UnPack,EvoId
nodes = 10
mysegment = DiscreteSegment(1,nodes)

# other possibility: defining the space as a graph
nodes = 10
g = LightGraphs.grid([nodes,1]) # cf LightGraphs.jl for options to generate a graph
mysegmentgraph = GraphSpace(g)
wholespace = (mysegmentgraph,)
using GraphPlot
# plotting the graph
gplot(g, collect(1:nodes), collect(1:nodes))

# Definition of birth and death rate
K0 = 1000 # We will have in total 1000 individuals
b(X,t) = 1 / nodes
d(X,Y,t) = (X[1] ≈ Y[1]) / K0
# Mutation / dispersal parameters
mu = [1.]
D = (1.5,)
# maximum size, tend
NMax = 2000
tend = 300.
# wrapping up all the parameters
p = Dict{String,Any}();@pack! p = D,mu,NMax

# definining world 0 and running
myagents = [Agent(wholespace,(5,),ancestors=true,rates=true) for i in 1:K0/nodes]
w0 = World(myagents,wholespace,p)
@time sim = run!(w0,Gillepsie(),tend,b,d)

### Plotting size of the world
myagents = [Agent(wholespace,(5,),ancestors=true,rates=true) for i in 1:K0/nodes]
w0 = World(myagents,wholespace,p) # we need to reinitialise the world
@time sim = run!(w0,Gillepsie(),tend,dt_saving=2.,b,d)
wsize = [length(w) for w in sim[:]]
using Plots
Plots.plot(get_tspan(sim),wsize,
                label = "",
                ylabel = "Metapopulation size",
                xlabel ="time",
                grid = false)
# savefig(joinpath(@__DIR__, "delta_comp_wsize.png"))

## plotting position through time
using Plots
Plots.plot(sim,
        label = "",
        ylabel = "Geographical position",
        grid = false,
        markersize = 10)
# savefig(joinpath(@__DIR__, "delta_comp_pos.png"))
