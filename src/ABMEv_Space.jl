using LightGraphs
"""
    abstract type AbstractSpace{Dim,T,F} end
`Dim` is the dimension of the space, `T` is the element type, `ife` is a bool which is `true`
when space is finite
"""

abstract type IsFinite{T} end

#ife stands for is finite
abstract type AbstractSpace{Dim,T,I} end
AbstractSpacesTuple = Tuple{Vararg{AbstractSpace}}
Base.ndims(x::AbstractSpace{Dim,T,I}) where {Dim,T,I} = Dim
Base.isfinite(x::AbstractSpace{Dim,T,IsFinite{t}}) where {Dim,T,t} = t #not sure we need this
Base.eltype(::AbstractSpace{Dim,T,I}) where {Dim,T,I} = T

SpaceType=Union{Nothing, AbstractSpace} # not sure what is this used for

struct GraphSpace{T} <: AbstractSpace{1,T,IsFinite{true}}
    g::AbstractGraph{T}
end

abstract type AbstractSegment{T}  <: AbstractSpace{1,T,IsFinite{true}} end
struct ContinuousSegment{T} <:  AbstractSegment{T}
    s::T
    e::T
end

struct DiscreteSegment{T} <: AbstractSegment{T}
    s::T
    e::T
end

struct Real1DSpace{T} <: AbstractSpace{1,T,IsFinite{false}} end

# TODO: find a way to put a type on get_inc
# TODO: there is probably a better way of dealing with get_inc

function get_inc(D,s::AbstractSpace{Dim,T,I}) where {Dim,T,I<:IsFinite{false}}
    if Dim > 1
        return D[:] .* rand(T,Dim)
    else
        return D * rand(T)
    end
end
get_inc(x,D,s::AbstractSpace{Dim,T,I}) where {Dim,T,I<:IsFinite{false}} = get_inc(D,s)

function get_inc(x,D,s::ContinuousSegment{T}) where {T}
    inc = D * rand(T)
    return reflect1D(x,inc,s)
end

function get_inc(x,D,s::DiscreteSegment{T}) where {T}
    inc = D * rand(T)
    return round(reflect1D(x,inc,s))
end

function get_inc(x,D,s::GraphSpace{T}) where {T}
    niter = round(rand(T) < D)
    inc = D * rand(T)
    return round(reflect1D(x,inc,s))
end

 """
     function increment_x!(a::Agent{MixedAgent,U},t::U,p::Dict) where U
 This function increments first trait of agent with integer values, that are automatically reflected between 1 and p["nodes"].
Other traits are incremented as usual.
TODO : make it work for a graph type landscape, where domain is not a line anymore.
 """
 function increment_x!(a::Agent{MixedAgent,U},t,p::Dict) where U
     tdim = length(p["D"])
     inc = [round(get_inc_reflected(get_x(a,1),p["D"][1] *randn(),1,p["nodes"]))]
     if  tdim > 1
         inc = vcat(inc,(rand(tdim-1) < p["mu"][2:end]) .* p["D"][2:end] .* randn(tdim-1))
     end
     a.x_history = hcat(a.x_history, get_x(a) + reshape(inc,:,1));
     push!(a.t_history,t)
end

"""
function get_inc_reflected(x::Number,inc::Number,s=-1,e=1)
    Here we increment the trajectory of trait 1 such that it follows a reflected brownian motion (1D)
"""
function reflect1D(x::Number,inc::Number,s<:AbstractSegment)
    if x + inc < s
        inc = 2 * ( s.s - x ) - inc
    elseif  x + inc > e
        inc = 2 * ( s.e - x ) - inc
    else
        return inc
    end
    reflect1D(x,inc,s)
end
