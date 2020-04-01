using ABMEv, Test, JLD2,Random

@testset "ABMEv" begin
    include("gillepsie.jl")
    include("wrightfisher.jl")
    include("metrics.jl")
end
