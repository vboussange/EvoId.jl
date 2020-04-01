cd(@__DIR__)
Random.seed!(0)
import ABMEv:update_rates_std!


a = 0;
sigma_K = .9;
sigma_a = .7;
K0 = 1000;
K(X) = gaussian(X[1],0.,sigma_K)
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
# α(X,Y) = 0.
p_default = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 1.5,
        "NMax" => Int(10000))
na_init = K0
world0 = new_world_G(na_init,p_default,spread = .01)
worldall,p_default["tspan"] = runWorld_store_G(p_default,world0)
world_alive_test = collect(skipmissing(worldall[:,end]))
# @save "gillepsie_test.jld2" world_alive
@load "gillepsie_test.jld2" world_alive
## Testing
@testset "Gillepsie Algorithm" begin
        @testset "Testing global functioning" begin
                @test size(worldall,2) > 1
                @test p_default["tspan"][end] >= p_default["tend"]
        end
        ## Comparing simulation
        xarray = get_xarray(world_alive,1);xarray_test = get_xarray(world_alive_test,1);
        @test xarray ≈ xarray_test

        @testset "Testing update rates matrix" begin
                bs_end = get_b.(world_alive);ds_end = get_d.(world_alive)
                update_rates_std!(world_alive,p_default,0.);
                bs_recalculated = get_b.(world_alive);ds_recalculated = get_d.(world_alive);

                @test bs_end ≈ bs_recalculated;

                @test ds_end ≈ ds_recalculated;

        end
end
