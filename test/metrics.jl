using LightGraphs
using Test
using Revise,ABMEv
using UnPack
myspace1 = (RealSpace{1,Float64}(),)
myspace2 = (RealSpace{1,Float64}(),RealSpace{1,Float64}())
myspace3 = (DiscreteSegment(Int16(1),Int16(10)),RealSpace{1,Float64}())
K0 = 1000; σ = 1e-1
a1 = [Agent(myspace1, (σ,)  .* randn() .- .5,ancestors=true) for i in 1:K0]
a2 = [Agent(myspace2,tuple((σ, σ)  .* randn(2) .- .5...),ancestors=true) for i in 1:K0]
a3 = [Agent(myspace3, (rand(Int16.(1:10)), 1e-1* randn() + 5.5 ),ancestors=true) for i in 1:K0]
D = (1.,);
mu = [1.,1.]
NMax = 1000
p1 = Dict{String,Any}();@pack! p1 = D,mu,NMax
D = (1.,1.);
p2 = Dict{String,Any}();@pack! p2 = D,mu,NMax
D = (Int16(0.),0.)
p3 = Dict{String,Any}();@pack! p3 = D,mu,NMax

w1 = World(a1,myspace1,p1)
w2 = World(a2,myspace2,p2)
w3 = World(a3,myspace3,p3)

## testing variance
@testset "Testing metrics" begin
    @testset "var" begin
        @test first(var(w1)) ≈ (σ).^2 atol=0.001
        @test first(var(w2,trait=2)) ≈ (σ).^2 atol=0.001
    end

    ## testing covgeo
    @testset "covgeo" begin
        # @test covgeo(w1) ≈ (σ).^2 atol=0.001
         for i in covgeo(w1,1)
             # @test i ≈ (σ).^2 atol=0.001
         end
     end
     # not sure this is the bestway of testing
     # there is a problem here
     @testset "covgeo2d" begin
         cmat = covgeo(w2,2);
         smat = [σ^2 0; 0 σ^2]
         # @test cmat ≈ smat atol=0.01
      end
      @testset "Alpha diversity" begin
          α = get_alpha_div(w3,2);
          @test abs(α) < Inf
       end
       @testset "Beta diversity" begin
           β = get_beta_div(w3,2);
           @test abs(β) < Inf
       end
       @testset "Isolation by history - hamming distance" begin
           a1 = Agent((DiscreteSegment(1,10),),[(1,),(2,),(3,)],[0,1,4],ancestors=true)
           a2 = Agent((DiscreteSegment(1,10),),[(1,),(10,),(3,),(10,)],[0,3,4,5],ancestors=true)
           @test get_dist_hist(a1,a2,(x,y)->y!=x,1) ≈ 3.0
       end
end

@testset "Geotrait computation" begin
    a = Agent(myspace3, (Int16.(1), randn() ),ancestors=true); increment_x!(a,myspace3,p3,2.0);
    @test get_geo(a,2.0) ≈ 2.0
end

# TODO needs to test hamming distance
