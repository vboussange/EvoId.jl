b(X,t) = gaussian(X[1],0.,sigma_K)
d(X,Y,t) = gaussian(X[1],Y[1],sigma_a)/K0
p = Dict{String,Any}();@pack! p = d,b
@testset "Utils" begin
    @test numargs(d) == 3
    @test numargs(b) == 2
end
