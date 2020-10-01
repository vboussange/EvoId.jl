# this used to be the worldalive

# TODO: do a constructor that ensures the parameters numerics are of the same type as the agents
mutable struct World{A<:AbstractAgent, S<:AbstractSpacesTuple,T<:Number}
    agents::Vector{AbstractAgentM}
    space::S
    parameters::Dict
    t::T
end

#constructor
function World(w::Vector{A},s::S,p::Dict,t::T=0.) where {A<:AbstractAgent,S<:AbstractSpacesTuple,T}
    if typeof(p["D"]) != eltype(skipmissing(w)[1])
        throw(ArgumentError("Diffusion coefficient does not match with underlying space\n `D::Tuple`"))
    end
    ww = vcat(w,repeat([missing],Int(p["NMax"] - length(w))))
    World{A,S,T}(ww,s,p,t)
end

# this throws an iterators of agents in the world
agents(world::World) = skipmissing(world.agents)
parameters(world::World) = world.parameters
time(w::World) = w.t
space(w::World) = w.space
maxsize(w::World) = w.parameters["NMax"]
_findfreeidx(w::World) = findfirst(ismissing,w.agents)
# this throws indices that are occupied by agents
_get_idx(world::World) = collect(eachindex(agents(world)))
# this throws agents of an abstract array of size size(world)
import Base:size,getindex
Base.size(world::World) = length(world.agents) - count(ismissing,world.agents)
Base.copy(w::W) where {W<:World} = W(copy.(w.agents),w.space,w.parameters,copy(w.t))
## Accessors
"""
$(SIGNATURES)
Get x of world without geotrait.
"""
Base.getindex(w::World,i::Int) = w.agents[_get_idx(w)[i]]

addAgent!(w::World,a::AbstractAgent) = begin
    idx = _findfreeidx(w)
    w.agents[idx] = a
    return nothing
end
removeAgent!(w::World,i::Int) = begin
    w.agents[_get_idx(w)[i]] = missing
    return nothing
end

update_clock!(w::World{A,S,T},dt) where {A,S,T} = begin
    w.t = convert(T,sum(w.t + dt))
    return nothing
end

#TODO : code it
"""
$(SIGNATURES)
Returns trait of every agents of world in the form of an array which dimensions corresponds to the input.
If `trait = 0` , we return the geotrait.
"""
get_x(w::World,trait) = vcat(getindex.(agents(w),trait)...)

"""
$(SIGNATURES)
Returns every traits of every agents of `world` in the form **of a one dimensional array** (in contrast to `get_x`).
If `geotrait=true` the geotrait is also added to the set of trait, in the last line.
If you do not want to specify `t` (only useful for geotrait), it is also possible to use `get_xarray(world::Array{T,1}) where {T <: Agent}`.
"""
function get_xarray(world::World,geotrait::Bool=false)
    xarray = hcat(get_x.(agents.(world))...)
    if geotrait
        xarray = vcat( xarray, get_geo.(agents(world),world.t)')
    end
    return xarray
end
@deprecate get_xarray(world,geotrait=false) get_x(world,Colon())
