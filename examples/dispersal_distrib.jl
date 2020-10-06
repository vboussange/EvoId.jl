using LightGraphs
using Test
using Revise,ABMEv
mysegmentgraph = GraphSpace(grid([10,1]))
mysegment = DiscreteSegment(1,10)

myplot = Plots.plot(grid = false)
for D in [1e-3,1e-2,1e-1,1e0,1e1]
    graphinc = [get_inc(5,D,mysegmentgraph) for i in 1:10000]
    seginc = [get_inc(5,D,mysegment) for i in 1:10000]
    using StatsPlots
    StatsPlots.density!(graphinc,label="Dispersal distribution for graph segment")
    StatsPlots.density!(seginc,label="Dispersal distribution for regular segment",linestyle=:dash)
end
myplot
