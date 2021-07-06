## Testing parallel implementation
using Random;Random.seed!(0)
cd(@__DIR__)
using Distributed;addprocs(exeflags="--project")
using Test
@everywhere begin
        using EvoId
        sigma_a = 1.251;
        K0 = 1000;
        K(X) = 1 - 0.125 * sum(X.^2)
        α(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
end
p = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 10.)
na_init = K0
agents = [Agent( [1e-2]  .* randn(1) .- .5) for i in 1:K0]

# @test nprocs() == 5

# time should be less than 1 sec on the second run
# using BenchmarkTools
@time (worldall_test,p["tspan"]) = runWorld_store_WF(p,agents,reflected=false)

@test size(worldall_test,2) == 10

## load to  compare simulation
# using JLD2
# @save "wrightfisher_test.jld2" worldall p
# @load "wrightfisher_test.jld2" worldall
# xarray = get_x(worldall,1); xarray_test = get_x(worldall_test,1)


# @test xarray ≈ xarray_test
