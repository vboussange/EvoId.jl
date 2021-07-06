using LightGraphs
using Test
using Revise,EvoId

##### SPACES #####
mysegment = DiscreteSegment(1,10)
mygraph = GraphSpace(SimpleGraph(10,10))
real2d = RealSpace{2,Float64}()
myline = RealSpace{1,Float16}()
mydiscreteline = NaturalSpace{1,Int8}()
mycontinuoussegment = ContinuousSegment(-1.,1.)
myspace = (mysegment,mygraph,myline)
myspace2 = (mysegment,mycontinuoussegment,real2d)

@testset "Space properties" begin
    # checking essential properties of spaces
    @test isfinite(mysegment) ≈ true
    @test isfinite(mygraph) ≈ true
    @test isfinite(myline) ≈ false
    @test ndims(real2d) ≈ 2
    @test isfinite(mycontinuoussegment) ≈ true
    @test typeof(myspace) <: AbstractSpacesTuple
    @test eltype(myspace2) == Tuple{Int64,Float64,Float64}

    # increment on infinite spaces
    @test EvoId.get_inc(0.,myline) ≈ (0.)
    @test EvoId.get_inc(0.,mydiscreteline) ≈ (0.)
    @test !(EvoId.get_inc(1.,myline) ≈ 0.)
    @test !(get_inc(1,1.,myline) ≈ 0.)
    # @test !(get_inc(1,1.,mydiscreteline) ≈ 0.)


    @test typeof(EvoId.get_inc([1.,0.],real2d)) == Vector{Float64}
    @test typeof(get_inc([1.,0.],[1.,0.],real2d)) == Vector{Float64}
    @test typeof(EvoId.get_inc([1.],real2d)) == Vector{Float64}
    @test typeof(EvoId.get_inc(1.,real2d)) == Vector{Float64}
    # EvoId._get_inc([1.],real2d)
    # EvoId.initpos(myspace2)


    # increment on finite spaces
    # checking if reflection works
    @testset "Reflection" begin
        @test mysegment.s - eps() < get_inc(5.,100.,mysegment) + 5. < mysegment.e + eps()
        @test mycontinuoussegment.s < get_inc(0.,100.,mycontinuoussegment) < mycontinuoussegment.e
        mysegment2 = DiscreteSegment(-1,1)
        @test EvoId._reflect1D(0.,2.0,mysegment2) ≈ .0
        @test EvoId._reflect1D(0.,-2.0,mysegment2) ≈ .0
        @test EvoId._reflect1D(0.,4.0,mysegment2) ≈ .0
        @test EvoId._reflect1D(0.,1.1,mysegment2) ≈ 1 - .1
    end

    #checking if graph works
    @test prod([get_inc(1,10,mygraph) + 1 ∈ vertices(mygraph.g) for i in 1:30])
    @test prod([get_inc(1,nothing,mygraph) + 1 ∈ vertices(mygraph.g) for i in 1:30])
end

@testset "Agent properties" begin
    ##### AGENTS #######
    a1 = Agent(myspace)
    a2 = Agent(myspace,ancestors = true)
    a3 = Agent(myspace,[1,1,1.])
    a4 = Agent(myspace2,[1,1.,[1.,1]],rates=true)
    a5 = Agent(myspace2,ancestors=true)
    # increment test
    p_myspace = Dict("D"=>Union{Float64,Float16}[1e-10,1e-10,Float16(1e-10)],"mu" =>Union{Float64,Float16}[1.,1.,Float16(1)] )
    p_myspace2 = Dict("D"=>[1,1,[1,1]],"mu" =>[1,1,[1,1]])

    # basic test
    @test typeof(a1) <: AbstractAgent
    @test eltype(a1) == eltype(myspace)
    @test eltype(a5) == eltype(myspace2)
    @test typeof(a1) <: AbstractAgent
    @test typeof(a1) == typeof(a3)

    old_a1 = deepcopy(a1)
    @test prod((get_x(old_a1) .≈ get_x(increment_x!(a1,myspace,p_myspace,0.))))
    # @test !prod((get_x(old_a1) .== get_x(increment_x!(a1,myspace,p_myspace,0.))))
    @test nancestors(increment_x!(a2,myspace,p_myspace,2.)) > 1
    @test !isnothing(increment_x!(a4,myspace2,p_myspace2,2.))
    @test !isnothing(increment_x!(a5,myspace2,p_myspace2,2.))
    @test typeof(EvoId._get_xinc(a2,myspace,p_myspace,0.)) == typeof(get_x(a2))

    # testing if get_xinc does not modify the value of a4
    a4 = Agent(myspace2,[1,1.,[1.,1]],rates=true)
    _xinc = EvoId._get_xinc(a4,myspace2,p_myspace2,0.)
    @test get_x(a4) == [1,1.,[1.,1]]
end

# @show EvoId._get_xinc(a4,myspace2,p_myspace2,0.)
@testset "Multidimensional spaces" begin
    myspace = (RealSpace{10,Float64}(),)
    a1 = Agent(myspace)
    p_myspace = Dict("D"=>[1e0],"mu" =>[1.] )
    increment_x!(a1,myspace,p_myspace,0.)
    @test !(get_x(a1)[1] == ones(10))

    # testing independence
    myspace = (RealSpace{1000,Float64}(),)
    a1 = Agent(myspace)
    p_myspace = Dict("D"=>[1e0],"mu" =>[fill(.1,1000)] )
    increment_x!(a1,myspace,p_myspace,0.)
    @test isapprox(sum(get_x(a1)[1] .== ones(1000)),900,atol = 30)

    myspace = (RealSpace{1000,Float64}(),)
    a1 = Agent(myspace)
    p_myspace = Dict("D"=>[rand(1000)],"mu" =>[fill(.1,1000)] )
    increment_x!(a1,myspace,p_myspace,0.)
    @test isapprox(sum(get_x(a1)[1] .== ones(1000)),900,atol = 30)
end
