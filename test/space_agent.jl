using LightGraphs
using Test
using Revise,ABMEv

##### SPACES #####
mysegment = DiscreteSegment(1,10)
mygraph = GraphSpace(SimpleGraph(10,10))
myline = RealSpace{2,Float64}()
mycontinuoussegment = ContinuousSegment(-1,1)
myspace = (mysegment,mygraph,myline,mycontinuoussegment)
myspace2 = (mysegment,mygraph,myline)

@test isfinite(mysegment) ≈ true
@test isfinite(mygraph) ≈ true
@test isfinite(myline) ≈ false
@test ndims(myline) == 2
@test isfinite(mycontinuoussegment) ≈ true
@test typeof(myspace) <: AbstractSpacesTuple

# still does not work
@test get_inc([0.],myline) ≈ (0.,0.)
@test get_inc([0.,0.],myline) ≈ (0.,0.)
# checking if reflection works
@test mysegment.s < get_inc(5.,100.,mysegment) + 5. < mysegment.e
@test mycontinuoussegment.s < get_inc(0.,100.,mycontinuoussegment) < mycontinuoussegment.e
#checking if graph works
@test get_inc(1,10,mygraph) + 1 ∈ vertices(mygraph.g)


##### AGENTS #######
a1 = Agent(myspace)
a1bis = Agent(myspace2)
a2 = Agent(myspace,ancestors = true)
a3 = Agent(myspace,(1,1,1.))
a4 = Agent(myspace,(1,1,1.),rates=true)
typeof(a1) <: AbstractAgent

agents = [a1,a2,a3,a4]
agentsm = [a1,a2,a3,a4,missing]

AgentBasedModel(agents,myspace)
AgentBasedModel(agentsm,myspace)
