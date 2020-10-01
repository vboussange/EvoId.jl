myspace1 = (RealSpace{1,Float64}(),)
myspace2 = (RealSpace{2,Float64}(),)
myspace3 = (DiscreteSegment(Int16(1),Int16(10)),RealSpace{1,Float64}())
K0 = 1000; σ = 1e-1
w1 = [Agent(myspace1, (σ,)  .* randn() .- .5) for i in 1:K0]
w2 = [Agent(myspace2,(tuple((σ, σ)  .* randn(2) .- .5...),)) for i in 1:K0]
w3 = [Agent(myspace3, (rand(Int16.(1:10)), 1e-1* randn() + 5.5 )) for i in 1:K0]
p = Dict("mu" => [1. 1.],"D" => [0. 0.], "nodes" =>10 )
## testing variance
@testset "Testing metrics" begin
    @testset "var" begin
        @test first(var(w1)) ≈ (σ).^2 atol=0.001
        @test first(var(w2,trait=2)) ≈ (σ).^2 atol=0.001
    end

    ## testing covgeo
    @testset "covgeo" begin
        @test covgeo(w1) ≈ (σ).^2 atol=0.001
         for i in covgeo(w1,1)
             @test i ≈ (σ).^2 atol=0.001
         end
     end
     # not sure this is the bestway of testing
     # there is a problem here
     @testset "covgeo2d" begin
         cmat = covgeo(w2,2);
         smat = [σ^2 0; 0 σ^2]
         @test cmat ≈ smat atol=0.01
      end
      @testset "Alpha diversity" begin
          α = get_alpha_div(w3,1.0,2);
          @test abs(α) < Inf
       end
       @testset "Beta diversity" begin
           β = get_beta_div(w3,1.0,2);
           @test abs(β) < Inf
       end
       @testset "Isolation by history - hamming distance" begin
           a1 = Agent{StdAgent,Float64}([1 2 3],[0,1,4],0.,1.)
           a2 = Agent{StdAgent,Float64}([1 10 3 10],[0,3,4,5],0.,1.)
           @test get_dist_hist(a1,a2,(x,y)->y!=x,1) ≈ 3.0
       end
end

@testset "Geotrait computation" begin
    a = Agent{MixedAgent}( Float16[1, randn()] ); increment_x!(a,1.,p);
    @test get_geo(a,2.0) ≈ 2.0
end

# TODO needs to test hamming distance
