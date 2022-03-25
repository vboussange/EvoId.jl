using EvoId
using LightGraphs
using Plots

nodes = 7
g = star_graph(nodes)
landscape = GraphSpace(g) # spatial space
θ = [rand([-0.5, 0.5]) for i in 1:nodes] # optimal trait values
# X[1][] is the geographical position
# X[2][] corresponds to the adaptive traits
traitspace = RealSpace(1) # phenotypic space
evolspace = (landscape, traitspace) # space over which individuals are structured

const K = 100. # carrying capacity
b(X,t) = 1. - 0.5 * (θ[Int(X[1][])] - X[2][])^2 # birth function
d(X,Y,t) = (X[1] ≈ Y[1]) ? 1. / K : 0. # death function
alg = Gillepsie() # update rule
NMax = 2000 # maximum number of individuals
D = [nothing, 5e-2] # dispersion coefficient
mu = [1e-1, 1e-1] # mutation / migration rate


myagents = [Agent(evolspace, [rand(1:nodes,1), randn(1) * D[2]]) for i in 1:K] # array containing the founder individuals
# random position on the graph
# random position on the trait space centered around 0
w0 = World(myagents, evolspace, D, mu, NMax, 0.) # the initial world, defined at time 0.

tend = 500. # time horizon
t_saving_cb = collect(range(0., tend, length=200)) # time step where callback function is called
cb(w) = Dict("N" => length(w)) # callback function

println("Running simulation with callback function")
@time sim = run!(w0, alg, tend, b, d, cb=cb, t_saving_cb = t_saving_cb)
Plots.plot(sim.tspan, sim["N"])

println("Running simulation with `dt_saving`")
myagents = [Agent(evolspace, [rand(1:nodes,1),randn(1) * D[2]]) for i in 1:K]
w0 = World(myagents, evolspace, D, mu, NMax, 0.)
@time sim = run!(w0, Gillepsie(), tend, b, d, dt_saving = 2.0)
# Plots.plot(sim, trait=2)
