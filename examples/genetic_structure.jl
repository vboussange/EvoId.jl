using ABMEv,UnPack

##### Genotype space#####
dim_neutr = 1000
magicprop = 523728 / 32896
g = SimpleGraph{Int16}(dim_neutr,round(Int16,dim_neutr * magicprop))
initnode = argmax(eigenvector_centrality(g)) # This is the central node the we will use to instantiate the populations
myspace = (DiscreteSegment(Int8(1),Int8(nodes)),GraphSpace(g)) # union of vector spaces

K0 = 1000
mu = [1.,1.]
b(X,t) = 1 / nodes
d(X,Y,t) = (X[1] ≈ Y[1]) / K0
D = (5e-1,1.4825)

NMax = 2000
# tend = 1.5
tend = 3000
p_default = Dict{String,Any}();@pack! p_default = d,b,NMax,mu
myagents = [Agent(myspace,(rand(Int8(1):Int8(nodes)),initnode),ancestors=true,rates=true) for i in 1:round(K0/nodes)]
w0 = World(myagents,myspace,p_default,0.)

@time sim = run!(w0,Gillepsie(),tend,dt_saving=3.)

using GraphPlot
# This is to plot with consistency
locs_x, locs_y = spring_layout(g)
