using EVOID
using LightGraphs
using Plots

nodes = 7
g = star_graph(nodes)
landscape = GraphSpace(g)
θ = [rand([-1,1]) for i in 1:nodes]
traitspace = RealSpace(1)
evolspace = (landscape,traitspace)

K = 150
b(X,t) = max(1 - 0.5 * (θ[X[1]] - X[2])^2,0)
# b(X,t) =1.
d(X,Y,t) = (X[1] ≈ Y[1]) / K

D = [nothing,5e-2]
mu = [1f-1,1f-1]

tend = 200.
t_saving_cb = collect(range(0.,tend,length=300))
cb(w) = Dict("N" => size(w))

p = Dict("NMax" => 2000,
        "D" => D,
        "mu" => mu
        )


myagents = [Agent(evolspace,[rand(1:nodes),randn() * D[2]]) for i in 1:K]
w0 = World(myagents,evolspace,p)

println("Running simulation with callback")
@time sim = run!(w0,Gillepsie(),tend,b,d,cb=cb,t_saving_cb = t_saving_cb)
Plots.plot(sim.tspan, sim["N"])

println("Running simulation with `dt_saving`")
myagents = [Agent(evolspace,[rand(1:nodes),randn() * D[2]]) for i in 1:K]
w0 = World(myagents,evolspace,p)
@time sim = run!(w0,Gillepsie(),tend,b,d,dt_saving=2.0)
Plots.plot(sim,trait=2)
