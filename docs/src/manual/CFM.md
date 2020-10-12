# CFM algorithm

Champagnat Ferriere Méléard algorithm described in [Champagnat and Ferriere founding article](https://linkinghub.elsevier.com/retrieve/pii/S0040580905001632). This algorithm similar to Gillepsie algorithms, excepts that it runs much faster!

Indeed, at every time step, only the fitness of the individual picked at random is evaluated. Thus this algorithm is of order ``K`` times more efficient.

In order to use it, you need to feed to the dictionary parameters `p` a constant `Cbar<:Real` that is the upperbound of the maximum of the sum of the birth and death rates (cf article).

## An example on how to use it
```julia
using ABMEv,UnPack,Plots
myspace = (RealSpace{1,Float64}(),)
σ_b = .9;
σ_d = .7;
K0 = 1000
b(X,t) = 1.
d(X,Y,t) = gaussian(X[1],Y[1],σ_d)/K0 / gaussian(X[1],0.,σ_b)
Cbar = b([0],0.) + d([0],[0],0.)
D = (1e-2,)
mu = [.1]
NMax = 2000
tend = 1500
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax,Cbar
myagents = [Agent(myspace,(1e-2 * randn(),)) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
@time sim = run!(w0,CFM(),tend,dt_saving = 4)
```

!!! warning "Development"
    CFM gives an approximate time step. As of now, we do not manage to obtain qualitatively the same results as the Gillepsie algorithm.

```@autodocs
Modules = [ABMEv]
Pages   = ["algo/ABMEv_CFM.jl"]
```
