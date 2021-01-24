using LightGraphs
using Test
using Revise,ABMEv
using UnPack
myspace1 = (RealSpace{1,Float64}(),)
myspace2 = (RealSpace{1,Float64}(),RealSpace{1,Float64}())
myspace3 = (DiscreteSegment(Int16(1),Int16(10)),RealSpace{2,Float64}())
g = SimpleGraph(1000,4000)
myspace4 = (RealSpace{1,Float64}(),GraphSpace(g),)

K0 = 1000; σ = 1e-1
a1 = [Agent(myspace1, σ  .* randn(1) .- .5,ancestors=true) for i in 1:K0]
a2 = [Agent(myspace2,σ  .* randn(2) .- .5,ancestors=true) for i in 1:K0]
a3 = [Agent(myspace3, [rand(Int16.(1:10)), 1e-1.* randn(2) .+ 5.5 ],ancestors=true) for i in 1:K0]
a4 = [Agent(myspace4, [1.,rand(Int64(1):Int64(1000))],ancestors=true) for i in 1:K0]

D = (1.,);
mu = [1.,1.]
NMax = 1000
p1 = Dict{String,Any}();@pack! p1 = D,mu,NMax
D = (1.,1.);
p2 = Dict{String,Any}();@pack! p2 = D,mu,NMax
D = (Int16(0.),(0.,0.))
p3 = Dict{String,Any}();@pack! p3 = D,mu,NMax
D = (0.,0.)
p4 = Dict{String,Any}();@pack! p4 = D,mu,NMax

w1 = World(a1,myspace1,p1)
w2 = World(a2,myspace2,p2)
w3 = World(a3,myspace3,p3)
w4 = World(a4,myspace4,p4)

## testing variance
@testset "Testing metrics" begin
    @testset "var" begin
        @test first(var(w1)) ≈ (σ).^2 atol=0.001
        @test first(var(w2,trait=2)) ≈ (σ).^2 atol=0.001
        @test first(var(w4,trait=1)) < Inf
        # @test first(var(w4,trait=2)) < Inf
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
          α = get_alpha_div(w4,1);
          @test abs(α) < Inf
          α = get_alpha_div(w3,2,false)
          @test length(α) == 10
       end
       @testset "Abundance" begin
           @test isapprox(get_local_abundance(w3), K0/10,atol = 10)
           N = get_local_abundance(w3,false)
           @test length(N) == 10
        end
       @testset "Beta diversity" begin
           β = get_beta_div(w3,2);
           @test abs(β) < Inf
           β = get_beta_div(w4,1);
           @test abs(β) < Inf
       end
       @testset "Isolation by history - hamming distance" begin
           a1 = Agent((DiscreteSegment(1,10),),[(1,),(2,),(3,)],[0,1,4],ancestors=true)
           a2 = Agent((DiscreteSegment(1,10),),[(1,),(10,),(3,),(10,)],[0,3,4,5],ancestors=true)
           @test get_dist_hist(a1,a2,(x,y)->y!=x,1) ≈ 3.0
       end
end

@testset "Geotrait computation" begin
    a = Agent(myspace3, [Int16.(1), randn(2) ],ancestors=true); increment_x!(a,myspace3,p3,2.0);
    @test get_geo(a,2.0) ≈ 2.0
end

# TODO needs to test hamming distance

## testing for real space of dimension d>1
multispace = (DiscreteSegment{Int8}(1,9),RealSpace{3,Float64}(),)
K0 = 10000;
multia = [Agent(multispace, [rand(Int8(1):Int8(10)),randn(3) ],ancestors=true) for i in 1:K0]
D = [Int8(1),fill(1.,3)]
mu = [1]
NMax = 1000
multip = Dict{String,Any}();@pack! multip = D,mu,NMax
multiw = World(multia,multispace,multip)
a = multia[10]
inc = get_inc(1,D[1],multispace[2])
@test inc[1] != inc[2] # checking that we have independent increments in each direction
@test isapprox(mean(var(multiw,trait = 2)),1.,atol=1e-2)
@test prod([isapprox(i,0.,atol=1e-1) for i in mean(multiw,trait = 2)])
@test isapprox(get_alpha_div(multiw,2),1.,atol = 1e-2)
@test isapprox(get_beta_div(multiw,2),0.,atol = 1e-2)

## some test to write metrics functions
# world = multiw;trait=2
# g = groupby(a->a[1],agents(world))
# v = [var(World(subw,space(world),parameters(world)),trait=trait) for subw in values(g)]
# h = vcat(v...)
# mean(mean(h,dims=1))
# mean(h)
# mean([mean(var(Float64.(get_x(World(subw,space(world),parameters(world)),trait)),corrected=false)) for subw in values(g)])
#
# [var(Float64.(get_x(World(subw,space(world),parameters(world)),trait)),corrected=false) for subw in values(g)]
# subw = collect(values(g))[1]
# var(Float64.(get_x(World(subw,space(world),parameters(world)),trait)),corrected=false)
# var(World(subw,space(world),parameters(world)),trait=trait)
# mean(World(subw,space(world),parameters(world)),trait=1)
#
# m = [mean(World(subw,space(world),parameters(world)),trait=trait) for subw in values(g)]
# h=vcat(m...)
# var(h,dims=1,corrected=false)
