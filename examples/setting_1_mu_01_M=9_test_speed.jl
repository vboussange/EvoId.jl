using Random
cd(@__DIR__)
using Dates
using JLD2
using EvoId
using LightGraphs, UnPack
using Revise, BenchmarkTools

######################
##### param def ######
######################
Tf = Float32 #Type for traits
Ti = Int32 # Type for graph nodes
Rtype = Float32 # Type for rates
TimeType = Float32 # type for time steps

nodes = Ti(9)
g = star_graph(nodes)
dim_neutr = 300
neutralspace = RealSpace{dim_neutr,Tf}()
const K1 = Rtype(150)
# good definition
@inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / K1 : 0f0
NMax = 2000
# tend = Tf(500.)
tend = TimeType(500.)
D = [nothing, fill(Tf(5e-2), dim_neutr)]
b(X,t) = 1f0
t_saving_cb = collect(range(0, tend, length=300))
myspace = (GraphSpace(g),neutralspace)
mu = [Tf(0.1), fill(Tf(1e-1), dim_neutr)]
p_default = Dict{String,Any}();
myagents = [Agent(myspace, [[rand(Ti(1):Ti(nodes))], D[2].* randn(Tf,dim_neutr)], Rtype) for i in 1:nodes*Int(K1)]
w0 = World(myagents, myspace, D, mu, NMax, t = TimeType(0))

@time sim = run!(w0, Gillepsie(), tend, b, d);
# [ Info: simulation stopped at t=500.0001760438999, after 1318480 steps
#  61.532121 seconds (27.88 M allocations: 31.510 GiB, 6.09% gc time, 0.21% compilation time)

#######################
####### Plotting ######
#######################
if false
    t_saving_cb = collect(range(0., tend, length=100))
    cb(w) = Dict("N" => length(w), "betau" => get_beta_div(w, 2),)
    @time sim = run!(w0, Gillepsie(), tend, b, d, t_saving_cb=t_saving_cb,cb=cb);
    using PyPlot
    fig, ax = plt.subplots(1,2)
    ax[1].plot(get_tspan(sim),sim["N"])
    ax[2].plot(get_tspan(sim),sim["betau"])
    gcf()
end