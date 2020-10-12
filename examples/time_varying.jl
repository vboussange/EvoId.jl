using ABMEv,UnPack


ω = 2* π / 150 # angular frequency
optimal_trait(t) = sin(ω * t)
tend = 300
Plots.plot(1:tend,optimal_trait,label = "Optimal trait",xlabel = "time")
# savefig(joinpath(@__DIR__, "optimal_trait.png"))

myspace = (RealSpace{1,Float64}(),)
sigma_K = 1.;
K0 = 1000 # We will have in total 1000 individuals
b(X,t) = gaussian(X[1],optimal_trait(t),sigma_K)
d(X,Y,t) = 1/K0
D = (5e-2,)
mu = [1.]
NMax = 2000
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax
myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
@time sim = run!(w0,Gillepsie(),tend,dt_saving=3.)

using Plots
Plots.plot(sim, ylabel = "Adaptive trait")
savefig(joinpath(@__DIR__, "time_varying_pop.png"))
