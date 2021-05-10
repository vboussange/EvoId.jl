@testset "Utils" begin
    b(X,t) = gaussian(X[1],0.,sigma_K)
    d(X,Y,t) = gaussian(X[1],Y[1],sigma_a)/K0
    @test numargs(d) == 3
    @test numargs(b) == 2
end
