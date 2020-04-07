
K0 = 1000; σ = 1e-1
agents1 = [Agent( [σ]  .* randn(1) .- .5) for i in 1:K0]
agents2 = [Agent( [σ σ]  .* randn(2) .- .5) for i in 1:K0]

## testing variance
@testset "Testing metrics" begin
    @test first(var(agents1)) ≈ (σ).^2 atol=0.001
    @test first(var(agents2,trait=2)) ≈ (σ).^2 atol=0.001

    ## testing covgeo
    @test covgeo(agents1) ≈ (σ).^2 atol=0.001
    @testset "covgeo" begin
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
end

# TODO needs to test hamming distance
