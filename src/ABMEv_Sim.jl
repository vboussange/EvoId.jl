abstract type AbstractAlg end

# this used to be  world
mutable struct Simulation{A<:AbstractAgent, S<:AbstractSpacesTuple,T<:Number}
    agentarray::Vector{Vector{AbstractAgent}}
    space::S
    tspan::Vector{T}
    df_agg::Vector{Dict}
    p::Dict{String,Any}
end

# callbacks has to be of the form (names" => String[],"aggregates" => Function)

"""
$(SIGNATURES)
"""
function Simulation(w0::World{A,S,T};cb = nothing) where {A,S,T}
    tspan = zeros(1)
    #agentarray is of size 2 at the beginning
    isnothing(cb) ? df_agg = Dict[] : df_agg = [cb(w0)]
    Simulation{A,S,T}([copy.(agents(w0))],space(w0),tspan,df_agg,parameters(w0))
 end

get_tend(s::Simulation) = s.tspan[end]
get_size(s::Simulation) = length(s.tspan)
get_tspan(s::Simulation) = s.tspan
get_world(s::Simulation,i) = World(s.agentarray[i],s.space,s.p,s.tspan[i])
Base.getindex(s::Simulation,i) = s.agentarray[i]
import Base.lastindex
Base.lastindex(s::Simulation) = get_size(s)

Base.getindex(s::Simulation,measure::String) = [agg[measure] for agg in s.df_agg]


function Base.show(io::IO, s::Simulation{A,S,T}) where {A,S,T}
     println(io, "Simulation with agents of type", A)
 end
# TODO: define two functions with signatures
# function add_entry!(s::Simulation{A,S,T,F},w::World) where {A,S,T,F<:Function}
# function add_entry!(s::Simulation{A,S,T,F},w::World) where {A,S,T,F<:Nothing}

"""
$(SIGNATURES)
Add `w` with callbacks `s.cb` to `s` if provided
"""
function add_entry!(s::Simulation{A,S,T},w::World,cb) where {A,S,T}
    push!(s.agentarray,copy.(agents(w)))
    push!(s.tspan,w.t)
    if !isnothing(cb)
        push!(s.df_agg,cb(w))
    end
    return nothing
end

function add_entry_cb_only!(s::Simulation{A,S,T},w::World,cb) where {A,S,T}
    push!(s.agentarray,[Agent(w.space)]) # pushing NA agents
    push!(s.tspan,w.t)
    if !isnothing(cb)
        push!(s.df_agg,cb(w))
    end
    return nothing
end

#TODO : code it
function get_xnt(s::Simulation;trait = 1)
    return [getindex.(wa,trait) for wa in s.agentarray],[fill(t,size(s[j])) for (j,t) in enumerate(s.tspan)]
end
# get_x(agentarray::Array{T},t,trait::Integer) where {T <: AbstractAgent} = reshape(hcat(get_x.(agentarray,t,trait)),size(agentarray,1),size(agentarray,2))
# @deprecate get_x(agentarray::Array{T},t::Number,trait::Integer)


#TODO: code it
# function world2df(world::Array{T,1},geotrait=false) where {T <: Agent}
#     xx = get_xarray(world)
#     dfw = DataFrame(:f => get_fitness.(world))
#     for i in 1:size(xx,1)
#         dfw[Meta.parse("x$i")] = xx[i,:]
#     end
#     if geotrait
#         dfw[:g] = get_geo.(world)
#     end
#     return dfw
# end
#
# """
#     world2df(world::Array{T,1},t::Number,geotrait = false) where {T <: Agent}
# Converts the array of agent world to a datafram, where each column corresponds to a trait of the
# agent, and an extra column captures fitness.
# Each row corresponds to an agent
# """
# function world2df(world::Array{T,1},t::Number,geotrait = false) where {T <: Agent}
#     xx = get_xarray(world)
#     dfw = DataFrame(:f => get_fitness.(world))
#     for i in 1:size(xx,1)
#         dfw[Meta.parse("x$i")] = xx[i,:]
#     end
#     if geotrait
#         dfw[:g] = get_geo.(world,t)
#     end
#     return dfw
# end
