

abstract type IsFinite{T} end

#ife stands for is finite
"""
$(TYPEDEF)
`Dim` is the dimension of the space,
`T` is the element type,
`I` to indicate finiteness
"""
abstract type AbstractSpace{Dim,T,I} end
AbstractSpacesTuple = Tuple{Vararg{AbstractSpace}}
import Base:ndims,isfinite,eltype
Base.ndims(x::AbstractSpace{Dim,T,I}) where {Dim,T,I} = Dim
Base.isfinite(x::AbstractSpace{Dim,T,IsFinite{t}}) where {Dim,T,t} = t #not sure we need this
Base.eltype(::AbstractSpace{Dim,T,I}) where {Dim,T,I} = T
Base.ndims(ss::AbstractSpacesTuple) = length(ss)
Base.eltype(ss::AbstractSpacesTuple) where {Dim,T,I} = Tuple{eltype.(ss)...}

SpaceType=Union{Nothing, AbstractSpace} # not sure what is this used for

# Static spaces
abstract type AbstractStatSpace{Dim,T,I} <: AbstractSpace{Dim,T,I} end

"""
$(TYPEDEF)
"""
struct GraphSpace{T} <: AbstractStatSpace{1,T,IsFinite{true}}
    g::AbstractGraph{T}
end

abstract type AbstractSegment{T<:Number}  <: AbstractStatSpace{1,T,IsFinite{true}} end

"""
$(TYPEDEF)
"""
struct ContinuousSegment{T<:AbstractFloat} <:  AbstractSegment{T}
    s::T
    e::T
end

"""
$(TYPEDEF)
"""
struct DiscreteSegment{T<:Integer} <: AbstractSegment{T}
    s::T
    e::T
end

"""
$(TYPEDEF)
A real space with dimension N and type T
"""
struct RealSpace{N,T} <: AbstractStatSpace{N,T,IsFinite{false}} end
RealSpace(N) = RealSpace{N,Float64}()
"""
$(TYPEDEF)
A natural space with dimension N and type T
"""
struct NaturalSpace{N,T} <: AbstractStatSpace{N,T,IsFinite{false}} end

## Increments - specialised function
# TODO: find a way to put a type on D in get_inc

"""
$(SIGNATURES)
Returns increment corresponding to space `s`
"""
get_inc(x,D,s::AbstractStatSpace,t) = get_inc(x,D,s) # this is defined to skip representation of t for following specialised methods
get_inc(x,D,s::AbstractSpace{Dim,T,I}) where {Dim,T,I<:IsFinite{false}} = get_inc(D,s) # This is defined to skip representation of x for spaces which do not use reflections.

function get_inc(D,s::AbstractSpace{Dim,T,I}) where {Dim,T<:AbstractFloat,I<:IsFinite{false}}
    if Dim > 1
        return D .* randn(T,Dim)
    else
        return D * randn(T)
    end
end

function get_inc(D,s::AbstractSpace{Dim,T,I}) where {Dim,T<:Integer,I<:IsFinite{false}}
    if Dim > 1
        return round.(T,D .*randn(Float32,Dim))
    else
        return round(D * randn(Float32))
    end
end

#TODO: there is probably a better way of dealing with those two functions
function get_inc(x,D,s::ContinuousSegment{T}) where {T}
    inc = D * randn(T)
    return _reflect1D(x,inc,s)
end

function get_inc(x,D,s::DiscreteSegment{T}) where {T}
    inc = D * randn()
    return round(T,_reflect1D(x,inc,s))
end

# normal dispersal kernel that gets truncated
function get_inc(x,D::Number,s::GraphSpace{T}) where {T}
    niter = round(Int,abs(D*randn())) + 1
    # here we add +1 since randomwalk(s.g,x,niter) returns x
    if niter > 0
        return last(randomwalk(s.g,x,niter)) - x
    else
        return 0
    end
end
# short range dispersal kernel, jump to neighbour node
function get_inc(x,D::Nothing,s::GraphSpace{T}) where {T}
    return last(randomwalk(s.g,x,2)) - x
end

## Dynamic spaces
abstract type AbstractDynSpace{Dim,T<:Number} <: AbstractSpace{Dim,T,IsFinite{true}} end
"""
$(TYPEDEF)
A dynamic graph space.

# Example
`DynGraphSpace(g,f)`
Function `f(t)` takes as argument time, and returns the index of the graph to pick at time `t` from array `g`
"""
struct DynGraphSpace{T<:Number} <: AbstractDynSpace{1,T}
    g::Vector{AbstractGraph{T}}
    f #update function
end
# This is surely not elegant, but we have not found an other way to do it yet
function DynGraphSpace(g::Array{A},f) where A <: AbstractGraph
     DynGraphSpace{eltype(g[1])}(g,f)
 end

"""
$SIGNATURES
Returns the graph correseponding to `d::DynGraphSpace` at time `t`
"""
get_graph(d::DynGraphSpace,t) = d.g[d.f(t)]


## Increments - specialised functions for dynamic graphs
# normal dispersal kernel that gets truncated
function get_inc(x,D::Number,d::DynGraphSpace{T},t) where {T}
    niter = round(Int,abs(D*randn())) + 1
    # here we add +1 since randomwalk(s.g,x,niter) returns x
    if niter > 0
        return last(randomwalk(get_graph(d,t),x,niter)) - x
    else
        return 0
    end
end

# short range dispersal kernel, jump to neighbour node
function get_inc(x,D::Nothing,d::DynGraphSpace{T},t) where {T}
    return last(randomwalk(get_graph(d,t),x,2)) - x
end

"""
$(SIGNATURES)
Here we increment the trajectory of trait 1 such that it follows a reflected brownian motion (1D)
"""
function _reflect1D(x::Number,inc::Number,s::AbstractSegment)
    if x + inc < s.s
        inc = 2 * ( s.s - x ) - inc
    elseif  x + inc > s.e
        inc = 2 * ( s.e - x ) - inc
    else
        return inc
    end
    _reflect1D(x,inc,s)
end
