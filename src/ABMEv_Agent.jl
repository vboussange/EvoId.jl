abstract type Ancestors{T} end
abstract type Rates{T} end
hasancestors(::Type{Ancestors{T}}) where {T} = T #not sure we need it
hasrates(::Type{Rates{T}}) where {T} = T # not sure we need it

abstract type AbstractAgent{A<:Ancestors,R<:Rates} end # tc for time contingency, fit for fitness coeff
AbstractAgentM = Union{Missing,AbstractAgent}
mutable struct Agent{A<:Ancestors,R<:Rates,T,U,V} <: AbstractAgent{A,R}
    # history of traits for geotraits
    pos::Array{T,1}
    # birth time of ancestors
    t_history::Array{U,1}
    # death rate
    d::V
    #birth rate
    b::V
end

# infers position type and zeros
function initpos(s::S) where {S<:AbstractSpacesTuple}
    Eltype = eltype.(s)
    Dims = ndims.(s)
    pos = tuple()
    for i in 1:length(Eltype)
        if Dims[i] > 1
            pos = (pos...,tuple(zeros(Eltype[i],Dims[i])...))
        else
            pos = (pos...,zero(Eltype[i]))
        end
    end
    Tuple{Eltype...},pos
end

function Agent(s::S;ancestors=false,rates=false) where {S  <: AbstractSpacesTuple}
    T,pos = initpos(s)
    t = ancestors ? [Float64(.0)] : [nothing]
    U =  ancestors ? Float64 : Nothing
    d = rates ?  Float64(.0) : nothing
    b = d
    V = rates ?  Float64 : Nothing
    Agent{Ancestors{ancestors},Rates{rates},T,U,V}([pos],t,d,b)
end

# here pos is provided
function Agent(s::S, pos::P;ancestors=false,rates=false) where {P,S  <: AbstractSpacesTuple}
    T = Tuple{eltype.(s)...}
    if !(T <: P)
        throw(ArgumentError("Position provided does not match with underlying space"))
    end
    t = ancestors ? Float64(.0) : [nothing]
    U =  ancestors ? Float64 : Nothing
    d = rates ?  Float64(.0) : nothing
    b = d
    V = rates ?  Float64 : Nothing
    Agent{Ancestors{ancestors},Rates{rates},T,U,V}([pos],t,d,b)
end

# TODO : implement pretty print

import Base.copy
copy(a::Agent{A,R,T,U,V}) where {A,R,T,U,V} = Agent{A,R,T,U,V}(copy(a.pos),copy(a.t_history),copy(a.d),copy(a.b))
copy(m::Missing) = missing
copy(n::Nothing) = nothing

#####################
###Agent accessors###
#####################

"""
    get_x(a::Agent)
Returns trait i of the agent
"""
get_x(a::AbstractAgent) = a.x_history[end]

"""
    function get_geo(a::Agent{U,T},t::Number) where {U,T}
Returns geotrait of agent `a` at time `t`
"""
function get_geo(a::Agent{Ancestors{true},R,T,U,V},t::Number) where {R,T,U,V}
    tarray = vcat(a.t_history[2:end],convert(U,t))
    tarray .-= a.t_history
    return sum(get_xhist(a,1) .* tarray)
end
# This method can acces geotrait, while the second not
"""
    get_x(a::Agent,t::Number,i::Integer)
Returns trait `i` of the agent.
Geotrait corresponds to dimension `i=0`.
"""

get_x(a::AbstractAgent,t::Number,i::Integer) = i > 0 ? a.x_history[end][Int(i)] : get_geo(a,t)
get_x(a::AbstractAgent,i::Integer) = a.x_history[end][Int(i)]
"""
    get_t(a::Agent) = a.t_history[end]
Get time when agent born.
"""
get_t(a::AbstractAgent) = a.t_history[end]
get_xhist(a::Agent,i::Number) = a.x_history[Int(i),:]
get_xhist(a::Agent) = a.x_history
get_thist(a::Agent) = a.t_history
get_d(a::Agent) = a.d
get_b(a::Agent) = a.b
get_fitness(a::Agent) = a.b - a.d
get_dim(a::Agent) = size(a.x_history,1)
get_nancestors(a::Agent) = size(a.x_history,2)

#####################
###World accessors###
#####################
"""
    get_x(world::Array{T},trait::Integer) where {T <: Agent}
Get x of world without geotrait.
"""
get_x(world::Array{T},trait::Integer) where {T <: Agent} = trait > 0 ? reshape(hcat(get_x.(world,trait)),size(world,1),size(world,2)) : throw(ErrorException("Not the right method, need `t` as an argument"))
"""
    get_x(world::Array{T},t::Number,trait::Integer) where {T <: Agent}
Returns trait of every agents of world in the form of an array which dimensions corresponds to the input.
If `trait = 0` , we return the geotrait.
"""
get_x(world::Array{T},t::Number,trait::Integer) where {T <: Agent} = trait > 0 ? reshape(hcat(get_x.(world,trait)),size(world,1),size(world,2)) : reshape(hcat(get_geo.(world,t)),size(world,1),size(world,2))

"""
    function get_xarray(world::Array{T,1}) where {T <: Agent}
Returns every traits of every agents of world in the form of an array
"""
function get_xarray(world::Array{T,1}) where {T <: Agent}
    return hcat(get_x.(world)...)
end
"""
    function get_xarray(world::Array{T,1},t::Number,geotrait::Bool=false) where {T <: Agent}
Returns every traits of every agents of `world` in the form **of a one dimensional array** (in contrast to `get_x`).
If `geotrait=true` the geotrait is also added to the set of trait, in the last line.
If you do not want to specify `t` (only useful for geotrait), it is also possible to use `get_xarray(world::Array{T,1}) where {T <: Agent}`.
"""
function get_xarray(world::Array{T,1},t::Number,geotrait::Bool=false) where {T <: Agent}
    xarray = hcat(get_x.(world)...)
    if geotrait
        xarray = vcat( xarray, get_geo.(world,t)')
    end
    return xarray
end

## Modifiers
"""
    function increment_x!(a::Agent{StdAgent,U},t::U,p::Dict) where U
This function increments agent by random numbers specified in p
ONLY FOR CONTINUOUS DOMAINS
"""
function increment_x!(a::AbstractAgent{A,R},s<:AbstractSpacesTuple,p<:Dict{String,Any},t<:T) where {A<:Ancestors{true},R,T}
    @unpack D,mu = p
    _x = get_x(a)
    inc = similar(_x)
    for (i,s) in enumerate(ss)
        if rand() < mu[i]
            inc[i] = get_inc.(_x[i],D[i],s)
        end
    end
    push!(a.t_history,t)
    a.x_history = push!(a.x_history, _x + inc);
end

function increment_x!(a::AbstractAgent{A,R},s<:AbstractSpacesTuple,p<:Dict{String,Any},t<:T) where {A<:Ancestors{false},R,T}
    @unpack D,mu = p
    _x = get_x(a)
    inc = get_inc.(_x,D,s)
    a.t_history = [t]
    a.x_history = [_x + inc];
end

"""
    function tin(t::Number,a::Number,b::Number)
if t in [a,b) returns 1. else returns 0
"""

function tin(t::Number,a::Number,b::Number)
    return t>=a && t<b ? 1. : 0.
end

function split_move(t)
    return .0 + 1/100*(t-20.)*tin(t,20.,120.) + tin(t,120.,Inf64)
end

function split_merge_move(t)
    return .0 + 1/30*(t-10.)*tin(t,10.,40.) + tin(t,40.,70.) + (1- 1/30*(t-70.))*tin(t,70.,100.)
end
