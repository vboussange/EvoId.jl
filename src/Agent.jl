abstract type AbstractAgent{X} end # tc for time contingency, fit for fitness coeff

mutable struct AgentwithAncestors{X} <: AbstractAgent{X}
    # history of traits for geotraits
    x_history::Vector{Vector{Vector{X}}}
    # birth time of ancestors
    t_history::Vector{X}
    # death rate
    d::X
    #birth rate
    b::X
end

mutable struct Agent{X} <: AbstractAgent{X}
    # traits
    x::Vector{Vector{X}}
    # death rate
    d::X
    # birth rate
    b::X
end

# infers position type and zeros
function _initpos(s::S) where {S<:AbstractSpacesTuple}
    _nd = ndims.(s)
    T = eltype.(s)
    Tprom = eltype(promote_type(T...))
    pos = Vector{Tprom}[]
    pos = [ones(Tprom,_nd[i]) for i in 1:length(_nd)]
    Tprom,pos
end

# default initialiser
"""
    Agent(s; ancestors=false)
    Agent(s, pos; ancestors=false)

Returns an `Agent` living on the underlying space `s` 
with initial position `pos`. If `pos` not provided, 
initialises agent with 0 values everywhere

# Keyword arguments
is required for the Gillepsie algorithm, but not for CFM algorithm
* `ancestors`. Set `ancestors=true` when you want to store ancestors traits.
"""
function Agent(s::S; ancestors=false) where {S  <: AbstractSpacesTuple}
    T, pos = _initpos(s)
    # promoting b, d, t_history to position
    _zero = zero(T)
    if ancestors
        return AgentwithAncestors([pos], [_zero], _zero, _zero)
    else
        return Agent(pos, _zero, _zero)
    end
end

function Agent(s::S, pos_t::Vector, t::Vector{U}; ancestors=false) where {S <: AbstractSpacesTuple, U <: AbstractFloat}
    T = eltype.(s)
    Tprom = eltype(promote_type(T...))
    @assert (length(pos_t) == length(t)) "length of `pos` should match length of `t`"
    @assert (length(pos_t[1]) == length(s)) "number of traits does not match with number of spaces"
    @assert all(length.(pos_t[1]) .== ndims.(s)) "number of traits does not match with number of dimension spaces"

    # it should be clear that if ancestors is false, pos_t has only one element
    ancestors ? nothing : @assert length(pos_t) == 1
    # pos_t[1] is trait at t[1]
    # pos_t[1][1] are the traits on space s[1]
    # pos_t[1][1][1] is the trait value at dim 1
    pos2_t = convert(Vector{Vector{Vector{Tprom}}}, pos_t)
    _zero = zero(eltype(promote_type(T...)))
    if ancestors
        return AgentwithAncestors(pos2_t, t, _zero, _zero)
    else
        return Agent(pos2_t[1], _zero, _zero)
    end
end

Agent(s, pos; ancestors=false) = Agent(s, [pos], [0.], ancestors=ancestors)


import Base:copy,show
# deepcopy is necessary for agents
Base.deepcopy(a::AgentwithAncestors) = AgentwithAncestors(deepcopy(a.x_history),copy(a.t_history),copy(a.d),copy(a.b))
Base.deepcopy(a::Agent) = Agent(deepcopy(a.x),copy(a.d),copy(a.b))

# this function only copies the trait history and time (x,t), and set birth and death rates to 0.
# useful for birth events
copyxt(a::Agent{X}) where {X} = Agent(copy(a.x), zero(X), zero(X))
copyxt(a::AgentwithAncestors{X}) where {X} = AgentwithAncestors(copy(a.x_history),copy(a.t_history),zero(X),zero(X))

# this has to be overloaded for Base.copy(a::Agent) to work properly
# TODO: not sure this is required anymore
Base.copy(m::Missing) = missing
Base.copy(n::Nothing) = nothing

#####################
###Agent accessors###
#####################
# for agents with ancestors, we retrieve the last index correpsonding to own traits 
# ( other indices correspond to ancestors )
Base.getindex(a::AgentwithAncestors, i) = a.x_history[end][i]
Base.getindex(a::Agent, i) = a.x[i]

"""
    get_x(a::Agent)
Returns trait i of the agent.
"""
get_x(a::AgentwithAncestors) = a.x_history[end]
get_x(a::Agent) = a.x

"""
    get_geo(a)
Returns geotrait of agent `a` at time `t`.
"""
function get_geo(a::AgentwithAncestors{X},t) where {X}
    tarray = vcat(a.t_history[2:end], t)
    tarray .-= a.t_history
    return sum(get_xhist(a,1) .* tarray)
end

"""
    get_t(a)
Get time when agent born.
"""
get_t(a::AgentwithAncestors) = a.t_history[end]

# TODO: change this with getindex
get_xhist(a::AgentwithAncestors,i) = [a.x_history[t][i] for t in 1:length(a.x_history)]

get_xhist(a::AgentwithAncestors) = a.x_history

get_thist(a::AgentwithAncestors) = a.t_history

get_d(a::AbstractAgent) = a.d

get_b(a::AbstractAgent) = a.b

get_fitness(a::AbstractAgent) = a.b - a.d

# TODO : we can surely extract N in Agent{A,R,Tuple{Vararg{S,N}},U,V}
# Inspiration : where U <: Union{Missing,Agent{T}} where T
Base.length(a::AgentwithAncestors) = length(a.x_history[end])
Base.length(a::Agent) = length(a.x)

nancestors(a::AgentwithAncestors) = length(a.x_history)

import Base.zero

Base.zero(t::Tuple{Vararg{Union{Number,Tuple{Vararg{Number}}}}}) = [zero.(e) for e in t]

import Base.(+)
(+)(t1::Tuple{Vararg{T,N}},t2::Tuple{Vararg{T,N}}) where {T<:Number,N}= tuple([t1[i] + t2[i] for i in 1:length(t1)]...)

function _get_xinc(a, s, D, mu, t)
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
    increment_x!(a, s, D, mu, t)
This function increments (i.e. mutates) agent `a` according to `D`, `mu`.
Note: for `a::AgentwithAncestors`, an extra row is added, since `increment_x!` is
called for birth events after copying the mum traits. Mum traits are therefore Ancestors traits.
"""
function increment_x!(a::AgentwithAncestors, s, D, mu, t)
    push!(a.t_history, t)
    push!(a.x_history, _get_xinc(a, s, D, mu, t))
    return a
end

function increment_x!(a::Agent, s, D, mu, t)
    a.x = _get_xinc(a, s, D, mu, t)
    return a
end

#################
###### IO #######
#################
function Base.show(io::IO, a::AbstractAgent{X}) where {X}
    println(io, "Agent with traits ", get_x(a))
end

Base.summary(A::AbstractAgent{X}) where {X} = string(TYPE_COLOR, nameof(typeof(a)), NO_COLOR," with trait type ",TYPE_COLOR, X)
