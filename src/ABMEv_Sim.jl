abstract type AbstractAlg end

# this used to be  world
mutable struct Simulation{A<:AbstractAgentM, S<:AbstractSpacesTuple,T<:Number,F}
    agentarray::Vector{AbstractAgentM}
    space::S
    tspan::Vector{T}
    cb::F
    df_agg
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
    !isnothing(cb.agg) ? df_agg = [Dict(cb.names .=> cb.agg(world))] : df_agg = nothing
    Simulation{A,S,T,cb.agg}([agentarray],space(w0),tspan,cb.agg,df_agg,parameters(w0))
 end

get_tend(s::Simulation) = s.tspan[end]
get_size(s::Simulation) = length(s.tspan)

"""
$(SIGNATURES)
Add `w` with callbacks `s.cb` to `s`
"""
function add_entry!(s::Simulation{A,S,T,F},w::World) where {A,S,T,F<:Function}
    push!(s.agentarray,copy.(collect(w.agents)))
    push!(s.tspan,w.t)
    push!(s.df_agg,Dict(s.cb.names .=> s.cb.agg(world)))
end

"""
$(SIGNATURES)
Add `w` to `s`
"""
function add_entry!(s::Simulation{A,S,T,F},w::World) where {A,S,T,F<:Nothing}
    push!(s.agentarray,copy.(collect(w.agents)))
    push!(s.tspan,w.t)
end

function world2df(world::Array{T,1},geotrait=false) where {T <: Agent}
    xx = get_xarray(world)
    dfw = DataFrame(:f => get_fitness.(world))
    for i in 1:size(xx,1)
        dfw[Meta.parse("x$i")] = xx[i,:]
    end
    if geotrait
        dfw[:g] = get_geo.(world)
    end
    return dfw
end

"""
    world2df(world::Array{T,1},t::Number,geotrait = false) where {T <: Agent}
Converts the array of agent world to a datafram, where each column corresponds to a trait of the
agent, and an extra column captures fitness.
Each row corresponds to an agent
"""
function world2df(world::Array{T,1},t::Number,geotrait = false) where {T <: Agent}
    xx = get_xarray(world)
    dfw = DataFrame(:f => get_fitness.(world))
    for i in 1:size(xx,1)
        dfw[Meta.parse("x$i")] = xx[i,:]
    end
    if geotrait
        dfw[:g] = get_geo.(world,t)
    end
    return dfw
end
