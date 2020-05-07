using Random;Random.seed!(0)
cd(@__DIR__)

## 1D Simulation
a = 0;
sigma_K = .9;
sigma_a = 1.251;
K0 = 1000;
# K(X) = gaussian(X[1],0.,sigma_K)
K(X) = 1 - 0.125 * sum(X.^2)
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
# α(X,Y) = 0.
p = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 10.)
na_init = K0
agents = [Agent( [1e-2]  .* randn(1) .- .5) for i in 1:K0]
worldall_test,p["tspan"] = runWorld_store_WF(p,agents);

## load to  compare simulation
# @save "wrightfisher_test.jld2" worldall p
@load "wrightfisher_test.jld2" worldall
xarray = get_x(worldall,1); xarray_test = get_x(worldall_test,1)
@testset "Wright Fisher Algorithm" begin
        @test xarray ≈ xarray_test
end
