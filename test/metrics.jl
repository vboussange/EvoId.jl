
K0 = 1000; σ = 1e-1
agents1 = [Agent( [σ]  .* randn(1) .- .5) for i in 1:K0]
agents2 = [Agent( [σ σ]  .* randn(2) .- .5) for i in 1:K0]
agentsd = [Agent{MixedAgent}( Float16[rand(1:10), 1e-1* randn() + 5.5] ) for i in 1:K0]
p = Dict("mu" => [1. 1.],"D" => [0. 0.], "nodes" =>10 )
## testing variance
@testset "Testing metrics" begin
    @testset "var" begin
        @test first(var(agents1)) ≈ (σ).^2 atol=0.001
        @test first(var(agents2,trait=2)) ≈ (σ).^2 atol=0.001
    end

    ## testing covgeo
    @testset "covgeo" begin
        @test covgeo(agents1) ≈ (σ).^2 atol=0.001
         for i in covgeo(agents1,1)
             @test i ≈ (σ).^2 atol=0.001
         end
     end
     # not sure this is the bestway of testing
     # there is a problem here
     @testset "covgeo2d" begin
         cmat = covgeo(agents2,2);
         smat = [σ^2 0; 0 σ^2]
         @test cmat ≈ smat atol=0.01
      end
      @testset "Alpha diversity" begin
          α = get_alpha_div(agentsd,1.0,2);
          @test abs(α) < Inf
       end
       @testset "Beta diversity" begin
           β = get_beta_div(agentsd,1.0,2);
           @test abs(β) < Inf
       end

end

@testset "Geotrait computation" begin
    a = Agent{MixedAgent}( Float16[1, randn()] ); increment_x!(a,1.,p);
    @test get_geo(a,2.0) ≈ 2.0
end

# TODO needs to test hamming distance
