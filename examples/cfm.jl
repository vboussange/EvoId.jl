using Revise,ABMEv,UnPack

myspace = (RealSpace{1,Float16}(),)
sigma_K = .9;
sigma_a = .7;
K0 = 5000
b(X,t) = gaussian(X[1],0.,sigma_K)
d(X,Y,t) = gaussian(X[1],Y[1],sigma_a)/K0
D = (Float16(1e-1),)
mu = [1.]
NMax = 7000
tend = 50
Cbar= b([0],0.) + d([0],[0],0.)
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax,Cbar
# myagents = [Agent(myspace,(- Float16(0.5) .+ .01 .* randn(Float16),)) for i in 1:K0]
myagents = [Agent(myspace,(0.,)) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
@time mysim = run!(w0,CFM(),tend,dt_saving=1.)

x,t = get_xnt(mysim,trait=1)
using Plots
Plots.scatter(t,x,color =:blue,labels="" )


########### GILLEPSIE
K0 = 1000
myspace = (RealSpace{1,Float64}(),)
sigma_K = .9;
sigma_a = .7;
b(X) = gaussian(X[1],0.,sigma_K)
d(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
D = (1e-1,)
mu = [1.]
NMax = 10000
tend = 100
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax #,Cbar
myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)

@time mysim = run!(w0,Gillepsie(),tend,dt_saving=1.)
x,t = get_xnt(sim,trait=1)
using Plots
Plots.scatter(t,x,color =:blue,labels="" )
