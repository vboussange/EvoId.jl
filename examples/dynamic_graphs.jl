# Dynamic graph

using UnPack,EvoId,LightGraphs
nodes = 10
g1 = LightGraphs.grid(Int8[9,1])
g2 = SimpleGraph(Int8(9))
g = [g1,g2] # array of graphs
# This is the function to implement DynGraphSpace.
# Note that it returns indices
function update_g(t)
    T = 100
    if sin(t/T*2*π) > 0
        1 # use graph g1
    else
        2 # use graph g2
    end
end
mydynamicgraph = DynGraphSpace(g,update_g)
wholespace = (mydynamicgraph,)

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
