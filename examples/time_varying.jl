###############################
# this example implements a birth rate that is time dependent
##############################

using EVOID,UnPack
using Plots

ω = 2* π / 150 # angular frequency
optimal_trait(t) = sin(ω * t)
tend = 300
Plots.plot(1:tend,optimal_trait,label = "Optimal trait",xlabel = "time")
# savefig(joinpath(@__DIR__, "optimal_trait.png"))

myspace = (RealSpace{1,Float64}(),)
sigma_K = 1.;
K0 = 300 # We will have in total 1000 individuals
b(X,t) = gaussian(X[1],optimal_trait(t),sigma_K)
d(X,Y,t) = 1/K0
D = (5e-2,)
mu = [1.]
NMax = 2000
p = Dict{String,Any}();@pack! p = D,mu,NMax
myagents = [Agent(myspace,(0,)) for i in 1:K0]
w0 = World(myagents,myspace,p)
@time sim = run!(w0,Gillepsie(),tend, b, d, dt_saving=3.)

Plots.plot(sim, ylabel = "Adaptive trait")
# savefig(joinpath(@__DIR__, "time_varying_pop.png"))
