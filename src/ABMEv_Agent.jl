abstract type Ancestors{T} end
abstract type Rates{T} end

abstract type AbstractAgent{A<:Ancestors,R<:Rates} end # tc for time contingency, fit for fitness coeff
AbstractAgentM = Union{Missing,AbstractAgent}
export AbstractAgentM

"""
$(TYPEDEF)
"""
mutable struct Agent{A<:Ancestors,R<:Rates,T<:Tuple,U,V} <: AbstractAgent{A,R}
    # history of traits for geotraits
    x_history::Array{T,1}
    # birth time of ancestors
    t_history::Array{U,1}
    # death rate
    d::V
    #birth rate
    b::V
end


eltype(a::Agent{A,R,T,U,V}) where {A,R,T,U,V} = T

# infers position type and zeros
function initpos(s::S) where {S<:AbstractSpacesTuple}
    Eltype = eltype.(s)
    Dims = ndims.(s)
    pos = tuple()
    for i in 1:length(Eltype)
        if Dims[i] > 1
            pos = (pos...,Eltype[i](ones(Dims[i])))
        else
            pos = (pos...,one(Eltype[i]))
        end
    end
    Tuple{Eltype...},pos
end

# default initialiser
"""
$(SIGNATURES)
    Initialises agent with 0 values everywhere
"""
function Agent(s::S;ancestors=false,rates=false) where {S  <: AbstractSpacesTuple}
    T,pos = initpos(s)
    t = 0.
    U =  Float64
    d = rates ?  Float64(.0) : nothing
    b = d
    V = rates ?  Float64 : Nothing
    Agent{Ancestors{ancestors},Rates{rates},T,U,V}([pos],[t],d,b)
end

# here pos is provided
"""
$(SIGNATURES)
    Initialises agent with `pos` provided
"""
function Agent(s::S, pos::P;ancestors=false,rates=false) where {P,S  <: AbstractSpacesTuple}
    T = eltype.(s)
    for (i,p) in enumerate(pos)
        if typeof(p) !== T[i]
            try
                p = convert(T[i],p)
            catch e
                throw(ArgumentError("Position provided does not match with underlying space"))
            end
        end
    end
    t = 0.
    U =  Float64
    d = rates ?  Float64(.0) : nothing
    b = d
    V = rates ?  Float64 : Nothing
    Agent{Ancestors{ancestors},Rates{rates},Tuple{T...},U,V}([pos],[t],d,b)
end

import Base:copy,show
Base.copy(a::A) where {A<:AbstractAgent} = A(copy(a.x_history),copy(a.t_history),copy(a.d),copy(a.b))

Base.copy(m::Missing) = missing

Base.copy(n::Nothing) = nothing

function Base.show(io::IO, a::AbstractAgent)
     print(io,"pos: ")
     show(io, a.x_history)
     println(io)
     print(io,"t: ")
     show(io, a.t_history)
 end

Base.summary(A::AbstractAgent) = string(TYPE_COLOR,nameof(typeof(a)),NO_COLOR," with uType ",TYPE_COLOR,eltype(a.x_history))

#####################
###Agent accessors###
#####################

"""
    get_x(a::Agent)
Returns trait i of the agent
"""

Base.getindex(a::Agent,i) = a.x_history[end][i]

get_x(a::Agent) = a.x_history[end]
@deprecate get_x(a) a[:]

"""
$(SIGNATURES)
Returns geotrait of agent `a` at time `t`
"""
function get_geo(a::Agent{A,R,T,U,V},t::Number) where {A<:Ancestors{true},R,T,U,V}
    tarray = vcat(a.t_history[2:end],convert(U,t))
    tarray .-= a.t_history
    return sum(get_xhist(a,1) .* tarray)
end

# This method can acces geotrait, while the second not
"""
$(SIGNATURES)
Returns trait `i` of the agent.
Geotrait corresponds to dimension `i=0`.
"""
get_x(a::Agent,t::Number,i::Integer) = i > 0 ? a.x_history[end][Int(i)] : get_geo(a,t)

"""
$(SIGNATURES)
Get time when agent born.
"""
get_t(a::Agent) = a.t_history[end]

# TODO: change this with getindex
get_xhist(a::AbstractAgent,i::Number) = [a.x_history[t][Int(i)] for t in 1:length(a.x_history)]

get_xhist(a::AbstractAgent) = a.x_history

get_thist(a::AbstractAgent) = a.t_history

get_d(a::AbstractAgent) = a.d

get_b(a::AbstractAgent) = a.b

get_fitness(a::AbstractAgent) = a.b - a.d

# TODO : we can surely extract N in Agent{A,R,Tuple{Vararg{S,N}},U,V}
# Inspiration : where U <: Union{Missing,Agent{T}} where T
Base.length(a::AbstractAgent) = length(a.x_history[end])

nancestors(a::Agent) = length(a.x_history)

import Base.zero

Base.zero(t::Tuple{Vararg{Union{Number,Tuple{Vararg{Number}}}}}) = [zero.(e) for e in t]

import Base.(+)
(+)(t1::Tuple{Vararg{T,N}},t2::Tuple{Vararg{T,N}}) where {T<:Number,N}= tuple([t1[i] + t2[i] for i in 1:length(t1)]...)

function _get_xinc(a::AbstractAgent,s::AbstractSpacesTuple,p::Dict,t::Number)
    @unpack D,mu = p
    _x = get_x(a)
    inc = zero(_x)
    for (i,ss) in enumerate(s)
        if rand() < mu[i]
            inc[i] = get_inc(_x[i],D[i],ss)
        end
    end
    tuple((_x .+ inc)...)
end

## Modifiers
"""
    $(SIGNATURES)
This function increments agent by random numbers specified in p
ONLY FOR CONTINUOUS DOMAINS
"""
function increment_x!(a::AbstractAgent{A,R},s::AbstractSpacesTuple,p::Dict,t::T) where {A<:Ancestors{true},R,T}
    push!(a.t_history,t)
    a.x_history = push!(a.x_history,_get_xinc(a,s,p,t))
    return a
end

function increment_x!(a::AbstractAgent{A,R},s::AbstractSpacesTuple,p::Dict,t::T) where {A<:Ancestors{false},R,T}
    a.t_history[1] = t
    a.x_history[1] = _get_xinc(a,s,p,t)
    return a
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
