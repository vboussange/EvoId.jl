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
    x_history::Vector
    # birth time of ancestors
    t_history::Vector{U}
    # death rate
    d::V
    #birth rate
    b::V
end

import Base:eltype
# This definition of eltype is to be discussed
eltype(a::Agent{A,R,T,U,V}) where {A,R,T,U,V} = T

# infers position type and zeros
function _initpos(s::S) where {S<:AbstractSpacesTuple}
    T = collect(eltype.(s))
    TT = collect(T)
    _nd = ndims.(s)
    for (i,n) in enumerate(_nd)
        if n > 1
            TT[i] = Vector{TT[i]}
        end
    end
    pos = Union{TT...}[]
    for i in 1:length(TT)
        if _nd[i] > 1
            push!(pos,ones(T[i],_nd[i]))
        else
            pos = push!(pos,one(T[i]))
        end
    end
    Tuple{T...},pos
end

# default initialiser
"""
$(SIGNATURES)
    Initialises agent with 0 values everywhere
    # args
    - `s` is the underlying space
    - `ancestors=true` when agents fitness needs to be updated at each time step.
    This is needed for the Gillepsie algorithm, but not for CFM algorithm
    - `ancestors=true` when one wants to store ancestors traits.
"""
function Agent(s::S;ancestors=false,rates=true) where {S  <: AbstractSpacesTuple}
    T,pos = _initpos(s)
    t = 0.
    U =  Float64
    d = rates ?  Float64(.0) : nothing
    b = d
    V = rates ?  Float64 : Nothing
    Agent{Ancestors{ancestors},Rates{rates},T,U,V}([pos],[t],d,b)
end

# here pos t should be provided
# this allows to initilise agents from any time steps knowing ancestors traits
function Agent(s::S,pos_t::Vector,t::Vector{U};ancestors=false,rates=true) where {S <: AbstractSpacesTuple, U <: AbstractFloat}
    T = eltype.(s)
    TT = collect(T) # we need an array to convert thereafter position
    _nd = ndims.(s)
    for (i,n) in enumerate(_nd)
        if n > 1
            TT[i] = Vector{T[i]}
        end
    end
    length(pos_t) == length(t) ? nothing : ArgumentError("length of `pos` should match length of `t`")
    # we convert all
    pos2_t = Vector{Union{TT...}}[]
    for pos in pos_t
        pos2 = Union{TT...}[]
        for (i,p) in enumerate(pos)
            try
                push!(pos2,convert(TT[i],p))
            catch e
                throw(ArgumentError("Position provided does not match with underlying space"))
            end
        end
        push!(pos2_t,pos2)
    end
    # U =  Float64
    d = rates ?  Float64(.0) : nothing
    b = d
    V = rates ?  Float64 : Nothing
    Agent{Ancestors{ancestors},Rates{rates},Tuple{T...},U,V}(pos2_t,t,d,b)
end

"""
Agent(s, pos;ancestors=false,rates=true)
    Initialises agent with initial position `pos` provided
    # args
    - `s` is the underlying space
    - `ancestors=true` when agents fitness needs to be updated at each time step.
    This is needed for the Gillepsie algorithm, but not for CFM algorithm
    - `ancestors=true` when one wants to store ancestors traits.
"""
Agent(s, pos; ancestors=false, rates=true) = Agent(s, [pos], [0.],ancestors=ancestors, rates=rates)


import Base:copy,show

Base.copy(a::A) where {A<:AbstractAgent} = A(copy(a.x_history),copy(a.t_history),copy(a.d),copy(a.b))

# this function only copies the trait history and time (x,t), and set birth and death rates to 0.
copyxt(a::Agent{A,R,T,U,V}) where {A,R,T,U,V<:Number} = Agent{A,R,T,U,V}(copy(a.x_history),copy(a.t_history),zero(V),zero(V))

copyxt(a::Agent{A,R,T,U,Nothing}) where {A,R,T,U} = Agent{A,R,T,U,Nothing}(copy(a.x_history),copy(a.t_history),nothing,nothing)

# this has to be overloaded for Base.copy(a::Agent) to work properly
Base.copy(m::Missing) = missing

Base.copy(n::Nothing) = nothing

function Base.show(io::IO, a::Agent{A,R,T,U,V}) where {A,R,T,U,V}
     println(io, "Agent with indices of type ", T)
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
get_x(a::AbstractAgent,t::Number,i::Integer) = i > 0 ? a.x_history[end][Int(i)] : get_geo(a,t)

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
    _x = deepcopy(get_x(a))
    for (i,ss) in enumerate(s)
        E = eltype(mu[i])
        S = eltype(ss)
        if length(mu[i]) > 1
            mut = (rand(E,ndims(ss)) .< mu[i]) .|> S
            _x[i] .+= mut .* get_inc(_x[i],D[i],ss,t)
        else
            mut = rand(E) < mu[i]
            if mut
                _x[i] += get_inc(_x[i],D[i],ss,t)
            end
        end
    end
    _x
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
