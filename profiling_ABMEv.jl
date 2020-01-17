using Distributed;addprocs()
@everywhere push!(LOAD_PATH,homedir()*"/ETHZ/projects/ABMEv.jl/src")
@everywhere using ABMEv,BenchmarkTools,SharedArrays

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
