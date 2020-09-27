using LightGraphs
mysegment = DiscreteSegment(1,10)
mygraph = GraphSpace(SimpleGraph(10,10))
myline = Real1DSpace{Float64}()

myspace = (mysegment,mygraph,myline)

@test isfinite(mysegment) ≈ true
@test isfinite(mygraph) ≈ true
@test isfinite(myline) ≈ false

a1 = Agent(myspace)
a2 = Agent(myspace,ancestors = true)
a3 = Agent(myspace,(1,1,1.))
a4 = Agent(myspace,(1,1,1.),rates=true)
typeof(a1) <: AbstractAgent

agents = [a1,a2,a3,a4]
agentsm = [a1,a2,a3,a4,missing]

AgentBasedModel(agents,myspace)
AgentBasedModel(agentsm,myspace)
