# using Test
# using Revise
# using ABMEv
# using UnPack
#
# b(X) = gaussian(X[1],0.,sigma_K)
# d(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
# p = Dict{String,Any}();@pack! p = d,b
# X = [.5]; Y = [.6];t = 0.
# _correct_timedep(p)
# @testset "Utils" begin
#     @test numargs(ABMEv.d) == 3
#     @test numargs(ABMEv.b) == 2
#     ABMEv.b(X,t) == b(X)
# end
