using Distributed;addprocs()
@everywhere push!(LOAD_PATH,homedir()*"/ETHZ/projects/ABMEv.jl/src")
@everywhere using ABMEv,BenchmarkTools,SharedArrays

## Testing update_afterbirth_std!
p= Dict("K0" => 1000.,
        "D" => [1e-2 - 1e-3],
        "mu" => [.1],
        "sigma_a" => [5e-1],
        "sigma_K" => [1e0],
        "n_alpha" => 2.,
        "n_K" => 2.,
        "tend" => 5.,
        "NMax" => Int(1e4))
na_init = 500
agent0 = [Agent([ .01  * randn() .- .5]) for i in 1:na_init]
world0 = vcat(agent0[:],repeat([missing],Int(p["NMax"] - na_init)))
C = SharedArray{Float64}((Int(p["K0"]),Int(p["K0"])))

@btime update_afterbirth_std!(skipmissing(world0),C,1,p)

## Testing get_inc_reflected
using ABMEv,BenchmarkTools
a = Agent(rand(1))
@btime get_inc_reflected(get_x(a)[1],.1 *randn())
@btime rand()
@btime get_x(a)[1]
@btime get_x(a,1)
@btime get_d(a)
@btime get_xhist(a)
@btime a.x_history[1,end]
@btime get_x(a)
@btime get_xhist(a)[:,end]

##
using ABMEv,BenchmarkTools
a = 0;
sigma_K = .9;
sigma_a = .7;
K(X) = gaussian(X[1],0.,sigma_K)
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/1000
# α(X,Y) = 0.
p_default = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 1000.,
        "NMax" => Int(10000))
na_init = 1000
world0 = new_world_G(na_init,p_default,spread = .01, offset = -.25)
tspan=zeros(1)
import ABMEv:update_rates_std!,updateWorld_G!
@btime update_rates_std!(skipmissing(world0),p_default,0.)
@btime updateWorld_G!(world0,p_default,update_rates_std!,tspan,true)
@btime update_afterbirth_std!(skipmissing(world0),1,p_default)
@btime update_afterdeath_std!(skipmissing(world0),[].8],p_default)
