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
tend = 100
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax

myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
w1 = copy(w0)
@info "Running simulation with Gillepsie algorithm"
@time sim = run!(w1,Gillepsie(),tend)
@test typeof(sim) <: Simulation

# @save "gillepsie_test.jld2" world_alive
# @load "gillepsie_test.jld2" world_alive
## Testing
@testset "Gillepsie Algorithm" begin
        @testset "Testing global functioning" begin
                @test size(sim,2) > 1
                @test get_tend(sim) >= tend
        end
        ## Comparing simulation
        # @testset "Matching new vs old results " begin
        #         xarray = get_x(w1,1);xarray_test = get_x(world_alive,1);
        #         @test prod(xarray .≈ xarray_test)
        # end

        @testset "Testing update rates matrix" begin
                @testset "birth event" begin
                        myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
                        w1 = World(myagents,myspace,p,0.);update_rates!(w1,Gillepsie())
                        mum_idx = 1
                        updateBirthEvent!(w0,Gillepsie(),1)
                        bs_end = get_b.(agents(w1));ds_end = get_d.(agents(w1))
                        update_rates!(w1,Gillepsie());
                        bs_recalculated = get_b.(agents(w1));ds_recalculated = get_d.(agents(w1))
                        @test prod(bs_end .≈ bs_recalculated)
                        @test prod(ds_end .≈ ds_recalculated)

                end

                @testset "death event" begin
                        myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
                        w1 = World(myagents,myspace,p,0.);update_rates!(w1,Gillepsie())
                        mum_idx = 1
                        updateDeathEvent!(w0,Gillepsie(),1)
                        bs_end = get_b.(agents(w1));ds_end = get_d.(agents(w1))
                        update_rates!(w1,Gillepsie());
                        bs_recalculated = get_b.(agents(w1));ds_recalculated = get_d.(agents(w1))
                        @test prod(bs_end .≈ bs_recalculated)
                        @test prod(ds_end .≈ ds_recalculated)
                end
        end
end
