using Revise,ABMEv,Test
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
world_alive = collect(skipmissing(world0))
@testset "testing global functioning" begin
        @test size(worldall,2) > 1
        @test p_default["tspan"][end] >= p_default["tend"]
end

@testset "testing update rates matrix" begin
        bs_end = get_b.(world_alive);ds_end = get_d.(world_alive)
        update_rates_std!(world_alive,p_default,0.);
        bs_recalculated = get_b.(world_alive);ds_recalculated = get_d.(world_alive);
                @testset "birth coefficient" begin
                        for i in 1:length(bs_end)
                        @test bs_end[i] ≈ bs_recalculated[i];
                end
                end
                @testset "death coefficient" begin
                        for i in 1:length(bs_end)
                        @test ds_end[i] ≈ ds_recalculated[i];
                end
                end
end
