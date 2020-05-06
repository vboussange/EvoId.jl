using ABMEv, Test, JLD2,Random

@testset "ABMEv" begin
    include("gillepsie.jl")
    # include("wrightfisher.jl")
    include("metrics.jl")
    # we might want to put this in a separate file at some point
    @testset "Reflection" begin
        @test get_inc_reflected(0.,2.0) ≈ .0
        @test get_inc_reflected(0.,-2.0) ≈ .0
        @test get_inc_reflected(0.,4.0) ≈ .0
        @test get_inc_reflected(0.,1.1) ≈ 1 - .1
    end
end
