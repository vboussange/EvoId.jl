abstract type AbstractAlg end

# this used to be  world
mutable struct Simulation{A<:AbstractAgent, S<:AbstractSpacesTuple,T<:Number,F}
    agentarray::Array{AbstractAgentM,2}
    space::S
    tspan::Vector{T}
    cb::NamedTuple
    df_agg::Vector{Dict}
    p::Dict{String,Any}
end

# callbacks has to be of the form (names" => String[],"aggregates" => Function)

"""
$(SIGNATURES)
"""
function Simulation(w0::World{A,S,T};cb=(names = String[],agg =nothing)) where {A,S,T}
    tspan = zeros(1);
    NMax = maxsize(w0)
    #agentarray is of size 2 at the beginning
    agentarray = copy.(collect(w0.agents))
    agentarray = hcat(agentarray,Array{Missing}(missing,NMax,1))
    !isnothing(cb.agg) ? df_agg = [Dict(cb.names .=> [f(w0) for f in cb.agg])] : df_agg = [Dict()]
    Simulation{A,S,T,typeof(cb.agg)}(agentarray,space(w0),tspan,cb,df_agg,parameters(w0))
 end

get_tend(s::Simulation) = s.tspan[end]
get_size(s::Simulation) = length(s.tspan)
get_tspan(s::Simulation) = s.tspan
Base.getindex(s::Simulation,i,j) = s.agentarray[i,j]
Base.getindex(s::Simulation,measure::String) = [agg[measure] for agg in s.df_agg]

import Base:lastindex,size
Base.lastindex(s::Simulation,i) = lastindex(s.agentarray,i)
Base.size(s::Simulation,i) = size(s.agentarray,i)

# TODO: define two functions with signatures
# function add_entry!(s::Simulation{A,S,T,F},w::World) where {A,S,T,F<:Function}
# function add_entry!(s::Simulation{A,S,T,F},w::World) where {A,S,T,F<:Nothing}

"""
$(SIGNATURES)
Add `w` with callbacks `s.cb` to `s` if provided
"""
function add_entry!(s::Simulation{A,S,T,F},w::World) where {A,S,T,F}
    i = get_size(s)
    j = size(s.agentarray,2)
    if i == j
        # we double the siwe of agent array
        s.agentarray = hcat(s.agentarray,Array{Missing}(missing,maxsize(w),j))
    end
    s.agentarray[:,i+1] .= copy.(collect(w.agents))
    push!(s.tspan,w.t)
    if !(F==Nothing)
        push!(s.df_agg,Dict(s.cb.names .=> [f(w) for f in s.cb.agg]))
    end
end

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
