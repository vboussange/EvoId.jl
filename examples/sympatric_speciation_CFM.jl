using EVOID,UnPack,Plots

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
dm = d([0],[0],0.);bm = 1.
p = Dict{String,Any}();@pack! p = dm,bm,D,mu,NMax
myagents = [Agent(myspace,(1e-2 * randn(Float64)),rates=false) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
@time sim = run!(w0,CFM(),tend,b,d,dt_saving = 10)


Plots.plot(sim,
        ylabel = "Adaptive trait",
        ylims = (-1,1),
        markersize = 2.)
