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
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax
myagents = [Agent(myspace,(1e-2 * randn(),),ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
@time sim = run!(w0,Gillepsie(),tend,dt_saving = 4)

using JLD2
@save joinpath(@__DIR__,"sim_sympatric_speciation.jld2") sim

Plots.plot(sim, ylabel = "Adaptive trait")
savefig(joinpath(@__DIR__, "sympatric_speciation.png"))

# plotting lineages
world = get_world(sim,get_size(sim))
xhistall = get_xhist.(world[:],1)
thist = get_thist.(world[:])
xplot = Plots.plot(thist,xhistall,
                linecolor = eth_grad_std[0.],
                label = "",
                # title = latexstring("\\sigma_\\mu=",@sprintf("%1.2f",world.p["D"][2][1]),", \\sigma_D=",@sprintf("%1.2f",world.p["D"][1])),
                grid = false,
                xlabel = "time",
                ylabel = "Historical adaptive trait"
                )a
savefig(joinpath(@__DIR__, "x_hist_sympatric_speciation.png"))
