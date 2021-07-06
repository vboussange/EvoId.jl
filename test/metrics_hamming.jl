cd(@__DIR__)
using Random
Random.seed!(0)
using LightGraphs
using Test
using Revise,EVOID
using UnPack,JLD2

myspace = (RealSpace{1,Float64}(),)
sigma_K = .9;
sigma_a = .7;
K0 = 1000;
b(X,t) = gaussian(X[1],0.,sigma_K)
d(X,Y,t) = gaussian(X[1],Y[1],sigma_a)/K0
D = [1e-2]
mu = [1.]
NMax = 10000
tend = 1.5
p = Dict{String,Any}();@pack! p = D,mu,NMax

myagents = [Agent(myspace,[0.,],ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
@info "Running simulation with Gillepsie algorithm"
@time sim = run!(w0,Gillepsie(),tend,b,d)

@testset "Hamming distances" begin
    @test typeof(get_xhist_mat(agents(w0))[1] )<: Array
    @test get_pairwise_average_isolation(w0) > 0
end
