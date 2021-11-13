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
Tf = Float32
Ti = Int32
nodes = Ti(9)
g = star_graph(nodes)
dim_neutr = 100
neutralspace = RealSpace{dim_neutr,Tf}()
K1 = 150
# good definition
@inbounds d(X, Y, t) = (X[1][] ≈ Y[1][]) ? 1f0 / Tf(K1) : 0f0
NMax = 2000
# tend = Tf(500.)
tend = Tf(200.)
D = [nothing, fill(Tf(5e-2), dim_neutr)]
b(X,t) = Tf(1.)
t_saving_cb = collect(range(0, tend, length=300))
myspace = (GraphSpace(g),neutralspace)
mu = [Tf(0.1),fill(Tf(1e-1),dim_neutr)]
p_default = Dict{String,Any}();
myagents = [Agent(myspace,Any[[rand(Ti(1):Ti(nodes))], D[2].* randn(Tf,dim_neutr)]) for i in 1:K1]
w0 = World(myagents, myspace, D, mu, NMax, t = Float32(0))

t_saving_cb = collect(range(0., tend, length=100))
cb(w) = Dict("N" => length(w), "betau" => get_beta_div(w, 2),)

@time sim = run!(w0, Gillepsie(), tend, b, d,t_saving_cb=t_saving_cb,cb=cb);

#######################
####### Plotting ######
#######################
using PyPlot
fig, ax = plt.subplots(1,2)
ax[1].plot(get_tspan(sim),sim["N"])
ax[2].plot(get_tspan(sim),sim["betau"])
gcf()