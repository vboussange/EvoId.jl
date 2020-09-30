cd(@__DIR__)
using Random
Random.seed!(0)
using LightGraphs
using Test
using Revise,ABMEv
using UnPack,JLD2

myspace = (RealSpace{1,Float64}(),)
sigma_K = .9;
sigma_a = .7;
K0 = 1000;
b(X) = gaussian(X[1],0.,sigma_K)
d(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
D = (1e-2,)
mu = [.1]
NMax = 10000
tend = 1.5
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax

myagents = [Agent(myspace,ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
w1 = copy(w0)

sim = run!(w0,Gillepsie(),tend)

world_alive_test = collect(skipmissing(sim[:,end]))
# @save "gillepsie_test.jld2" world_alive
@load "gillepsie_test.jld2" world_alive
## Testing
@testset "Gillepsie Algorithm" begin
        @testset "Testing global functioning" begin
                @test size(sim,2) > 1
                @test get_tend(sim) >= tend
        end
        ## Comparing simulation
        @testset "Matching new vs old results " begin
                xarray = get_x(world_alive,1);xarray_test = get_x(world_alive_test,1);
                @test xarray ≈ xarray_test
        end

        @testset "Testing update rates matrix" begin
                bs_end = get_b.(world_alive);ds_end = get_d.(world_alive)
                update_rates_std!(world_alive,p_default,0.);
                bs_recalculated = get_b.(world_alive);ds_recalculated = get_d.(world_alive);

                @test bs_end ≈ bs_recalculated;

                @test ds_end ≈ ds_recalculated;

        end
end
