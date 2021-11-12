using Random
length(ARGS) > 0 ? seed = (parse(Int,ARGS[1])) : seed = 1
Random.seed!(seed)
cd(@__DIR__)
using Dates
using JLD2
using EvoId
using LightGraphs, UnPack
using Revise, BenchmarkTools

######################
##### param def ######
######################
Tf = Float32
Ti = Int32
nodes = Ti(9)
g = star_graph(nodes)
dim_neutr = 300
neutralspace = RealSpace{dim_neutr,Tf}()
K1 = 150
# good definition
@inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / Tf(K1) : 0f0
# bad writing
# d(X, Y, t) = (X[1] ≈ Y[1]) / K1
NMax = 2000
# tend = Tf(500.)
tend = Tf(20.)
D = [nothing, fill(Tf(5e-2), dim_neutr)]
b(X,t) = Tf(1.)
t_saving_cb = collect(range(0, tend, length=300))
myspace = (GraphSpace(g),neutralspace)
mu = [Tf(0.1),fill(Tf(1e-1),dim_neutr)]
p_default = Dict{String,Any}();

myagents = [Agent(myspace,Any[[rand(Ti(1):Ti(nodes))], D[2].* randn(Tf,dim_neutr)]) for i in 1:nodes*K1]
w0 = World(myagents, myspace, D, mu, NMax)
@time s = run!(w0, Gillepsie(), tend, b, d);