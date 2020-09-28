##### SPACES #####
using LightGraphs
using Test
using Revise,ABMEv

mysegment = DiscreteSegment(1,10)
mygraph = GraphSpace(SimpleGraph(10,10))
myline = RealSpace{2,Float64}()
mycontinuoussegment = ContinuousSegment(-1,1)
myspace = (mysegment,mygraph,myline,mycontinuoussegment)
myspace2 = (mysegment,mygraph,myline)

@test isfinite(mysegment) ≈ true
@test isfinite(mygraph) ≈ true
@test isfinite(myline) ≈ false
@test isfinite(mycontinuoussegment) ≈ true
@test typeof(myspace) <: AbstractSpacesTuple



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
