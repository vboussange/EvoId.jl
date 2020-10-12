using Revise,ABMEv
nodes = 10
mysegment = DiscreteSegment(1,nodes)

using ABMEv, LightGraphs
nodes = 10
g = grid([nodes,1])
mysegmentgraph = GraphSpace(g)
wholespace = (mysegmentgraph,)
using GraphPlot
gplot(g, collect(1:nodes), collect(1:nodes))

K0 = 1000 # We will have in total 1000 individuals
b(X,t) = 1 / nodes
d(X,Y,t) = (X[1] ≈ Y[1]) / K0
mu = [1.]
D = (1.5,)

using UnPack
NMax = 2000
tend = 300.
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax
myagents = [Agent(wholespace,(5,),ancestors=true,rates=true) for i in 1:K0/nodes]
w0 = World(myagents,wholespace,p,0.)
@time sim = run!(w0,Gillepsie(),tend)

### Plotting size of the world
myagents = [Agent(wholespace,(5,),ancestors=true,rates=true) for i in 1:K0/nodes]
w0 = World(myagents,wholespace,p,0.) # we need to reinitialise the world
@time sim = run!(w0,Gillepsie(),tend,dt_saving=2.)
wsize = [length(w) for w in sim[:]]
using Plots
Plots.plot(get_tspan(sim),wsize,
                label = "",
                ylabel = "Metapopulation size",
                xlabel ="time",
                grid = false)
savefig(joinpath(@__DIR__, "delta_comp_wsize.png"))

## plotting position through time
using Plots
Plots.plot(sim,
        label = "",
        ylabel = "Geographical position",
        grid = false,
        markersize = 10)
savefig(joinpath(@__DIR__, "delta_comp_pos.png"))
