using LightGraphs
using Test
using Revise,ABMEv

##### SPACES #####
mysegment = DiscreteSegment(1,10)
mygraph = GraphSpace(SimpleGraph(10,10))
real2d = RealSpace{2,Float64}()
myline = RealSpace{1,Float16}()
mycontinuoussegment = ContinuousSegment(-1.,1.)
myspace = (mysegment,mygraph,myline)
myspace2 = (mysegment,mycontinuoussegment,real2d)

# checking essential properties of spaces
@test isfinite(mysegment) ≈ true
@test isfinite(mygraph) ≈ true
@test isfinite(myline) ≈ false
@test ndims(real2d) ≈ 2
@test isfinite(mycontinuoussegment) ≈ true
@test typeof(myspace) <: AbstractSpacesTuple
@test eltype(myspace2) == Tuple{Int64,Float64,Tuple{Float64,Float64}}

# increment on infinite spaces
@test get_inc(0.,myline) ≈ (0.)
@test !(get_inc(1.,myline) ≈ 0.)
@test !(get_inc(1,1.,myline) ≈ 0.)

@test typeof(get_inc([1.,0.],real2d)) == Tuple{Float64,Float64}
@test typeof(get_inc([1.,0.],[1.,0.],real2d)) == Tuple{Float64,Float64}
@test typeof(get_inc([1.],real2d)) == Tuple{Float64,Float64}
@test typeof(get_inc(1.,real2d)) == Tuple{Float64,Float64}
get_inc([1.],real2d)
ABMEv.initpos(myspace2)

# increment on finite spaces
# checking if reflection works
@test mysegment.s < get_inc(5.,100.,mysegment) + 5. < mysegment.e
@test mycontinuoussegment.s < get_inc(0.,100.,mycontinuoussegment) < mycontinuoussegment.e
#checking if graph works
@test prod([get_inc(1,10,mygraph) + 1 ∈ vertices(mygraph.g) for i in 1:30])


##### AGENTS #######
a1 = Agent(myspace)
a2 = Agent(myspace,ancestors = true)
a3 = Agent(myspace,(1,1,1.))
a4 = Agent(myspace2,(1,1,(1.,1)),rates=true)
a5 = Agent(myspace2,ancestors=true)

# basic test
@test typeof(a1) <: AbstractAgent
@test eltype(a1) == eltype(myspace)
@test eltype(a5) == eltype(myspace2)
@test typeof(a1) <: AbstractAgent

# increment test
p_myspace = Dict("D"=>[1,1,1],"mu" =>[1,1,1] )
p_myspace2 = Dict("D"=>[1,1,1],"mu" =>[1,1,1])
old_a1 = copy(a1)
@test !prod((get_x(old_a1) .≈ get_x(increment_x!(a1,myspace,p_myspace,0.))))
@test nancestors(increment_x!(a2,myspace,p_myspace,2.)) > 1
@test !isnothing(increment_x!(a4,myspace2,p_myspace2,2.))
@test !isnothing(increment_x!(a5,myspace2,p_myspace2,2.))
