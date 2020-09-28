# this used to be the worldalive
mutable struct World{A<:AbstractAgent, S<:AbstractSpacesTuple,T<:Number}
    agents::Vector{AbstractAgentM}
    space::S
    parameters::Dict{String,Any}
    t::T
end

parameters(world::World) = world.parameters
time(w::World) = w.time
space(w::World) = w.space
# this throws indices that are occupied by agents
_get_idx(world::World) = collect(eachindex(skipmissing(world.agents)))
# this throws agents of an abstract array of size size(world)
Base.getindex(w::World,i::Int) = w.agents[_get_idx(w)[i]]

# this throws an iterators of agents in the world
agents(world) = skipmissing(world.agents)
size(world) = size(world.agents) - count(ismissing,world.agents)
maxsize(w::World) = length(w.agents)
_findfreeidx(w::World) = findfirst(findfirst(ismissing,w.agents))
addAgent!(w::World,a::AbstractAgent) = begin
    idx = _findfreeidx(w)
    w.agents[idx] = a
    return nothing
end
removeAgent!(w::World,i::Int) = begin
    world[i] = missing
    return nothing
end

update_clock!(w::World{A,S,T},dt) where {A,S,T} = begin
    w.t = convert(T,sum(w.t + dt))
    return nothing
end
