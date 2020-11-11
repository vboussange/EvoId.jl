using ABMEv,UnPack,Plots

myspace = (RealSpace{1,Float64}(),)
σ_b = .9;
σ_d = .7;
K0 = 1000
b(X,t) = 1.
d(X,Y,t) = gaussian(X[1],Y[1],σ_d)/K0 / gaussian(X[1],0.,σ_b)
D = (1e-2,)
mu = [.1]
NMax = 2000
tend = 1500
p = Dict{String,Any}();@pack! p = D,mu,NMax
myagents = [Agent(myspace,(1e-2 * randn(Float64),),rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
@time sim = run!(w0,Gillepsie(),tend,dt_saving = 10,b,d)

using JLD2
@load joinpath(@__DIR__,"sim_sympatric_speciation.jld2") sim

Plots.plot(sim,
        ylabel = "Adaptive trait",
        ylims = (-1,1),
        markersize = 2.)
savefig(joinpath(@__DIR__, "sympatric_speciation.png"))
