
# this used to be  world
mutable struct Simulation{A<:AbstractAgentM, S<:AbstractSpacesTuple,T<:Number,N<:Number}
    worldarray::Array{AbstractAgentM,2}
    tspan::Vector{T}
    aggregates::Vector{Dict{String,Number}}
    p::Dict{String,Any}
end

function Simulation()
    tspan = zeros(1);
    worldall = reshape(copy.(world0),N,1);
    worldall = hcat(worldall,Array{Missing}(missing,N,1))

 end

get_tend(s::Simulation) = s.tspan[end]

function add_entry!(s::Simulation,w::World)
    j+=1;sw = size(worldall,2);
    # we use <= because we save at the end of the wile loop
    if sw <= j
        # we double the size of worldall
        worldall = hcat(worldall,Array{Missing}(missing,N,sw))
    end
    worldall[1:Int(N - count(ismissing,world0)),j] .= copy.(collect(skipmissing(world0)));
    push!(tspanarray,t)
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
