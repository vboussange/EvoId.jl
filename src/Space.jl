

abstract type IsFinite{T} end

# `Dim` is the dimension of the space,
# `T` is the element type,
# `I` to indicate finiteness

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
    GraphSpace(g)

Creates a Graph Space.

# Example
```julia
using LightGraphs
g = star_graph(7)
GraphSpace(g)
```

"""
struct GraphSpace{T} <: AbstractStatSpace{1,T,IsFinite{true}}
    g::AbstractGraph{T}
end

abstract type AbstractSegment{T<:Number}  <: AbstractStatSpace{1,T,IsFinite{true}} end

"""
    ContinuousSegment(s, e)

Creates a segment space, where individuals are reflected at both ends.

# Arguments
* `s` start of the segment
* `e` end of the segment
# Example
```julia
ContinuousSegment(1., 2.)
```
"""
struct ContinuousSegment{T<:AbstractFloat} <:  AbstractSegment{T}
    s::T
    e::T
end

"""
    DiscreteSegment(s, e)

Creates a discrete segement space, where individuals are reflected at both ends.

# Arguments
* `s` start of the segment
* `e` end of the segment

# Example

```julia
    DiscreteSegment(1, 2)
```
"""
struct DiscreteSegment{T<:Integer} <: AbstractSegment{T}
    s::T
    e::T
end

"""
    RealSpace{N,T}()

Creates a real space.

# Arguments
* `N` dimension of the space 
* `T` type of the underlying traits.
"""
struct RealSpace{N,T} <: AbstractStatSpace{N,T,IsFinite{false}} end
RealSpace(N) = RealSpace{N,Float64}()

"""
    NaturalSpace{N,T}
A natural space with dimension `N` and type `T`
"""
struct NaturalSpace{N,T} <: AbstractStatSpace{N,T,IsFinite{false}} end

## Increments - specialised function
# TODO: find a way to put a type on D in get_inc

"""
    get_inc(x, D, s)

Returns increment corresponding to space `s`
"""
get_inc(x, D, s::AbstractStatSpace,t) = get_inc(x,D,s) # this is defined to skip representation of t for following specialised methods
get_inc(x, D, s::AbstractSpace{Dim,T,I}) where {Dim,T,I<:IsFinite{false}} = get_inc(D,s) # This is defined to skip representation of x for spaces which do not use reflections.

function get_inc(D::DType, s::AbstractSpace{Dim,T,I}) where {DType, Dim,T<:AbstractFloat,I<:IsFinite{false}}
    return D .* randn(eltype(DType),Dim)
end


function get_inc(D::DType, s::AbstractSpace{Dim,T,I}) where {DType, Dim, T<:Integer, I<:IsFinite{false}}
    return round.(D .*randn(eltype(DType),Dim))
end

#TODO: there is probably a better way of dealing with those two functions
function get_inc(x,D::DType,s::ContinuousSegment{T}) where {T, DType}
    inc = D .* randn(T,1)
    return _reflect1D(x,inc,s)
end

function get_inc(x,D::DType,s::DiscreteSegment{T}) where {T,DType}
    inc = D .* randn(eltype(DType),1)
    return round.(_reflect1D(x,inc,s))
end

function get_inc(x,D::Nothing,s::DiscreteSegment{T}) where {T}
    inc = rand([one(T),-one(T)])
    return round.(T,_reflect1D(x,inc,s))
end

# normal dispersal kernel that gets truncated
function get_inc(x::xType, D::DType, s::GraphSpace{T}) where {xType, T, DType}
    niter = round(T,abs(D[1] * randn(eltype(DType)))) + 1
    # here we add +1 since randomwalk(s.g,x,niter) returns x
    if niter > 0
        inc = [last(randomwalk(s.g, convert(T,x[]), niter))] - x
        return  inc
    else
        return zero(x)
    end
end
# short range dispersal kernel, jump to neighbour node
function get_inc(x::xType, D::Nothing, s::GraphSpace{T}) where {xType,T}
    inc = [last(randomwalk(s.g, convert(T,x[]), 2))] - x
    return  convert(xType,inc) 
end

## Dynamic spaces
abstract type AbstractDynSpace{Dim,T<:Number} <: AbstractSpace{Dim,T,IsFinite{true}} end
"""
    DynGraphSpace(g, f)

A dynamic graph space.

# Arguments
* `g` the underlying graph
* `f` a function that takes as argument time,
and returns the index of the graph to pick at time `t` from array `g`

# Example
`DynGraphSpace(g,f)`
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
    get_graph(d, t)

Returns the graph correseponding to `d::DynGraphSpace` at time `t`
"""
get_graph(d::DynGraphSpace,t) = d.g[d.f(t)]


## Increments - specialised functions for dynamic graphs
# normal dispersal kernel that gets truncated
function get_inc(x, D::DType, d::DynGraphSpace{T},t) where {DType, T}
    niter = round(T,abs(D*randn(eltype(DType)))) + 1
    # here we add +1 since randomwalk(s.g,x,niter) returns x
    if niter > 0
        inc = [last(randomwalk(get_graph(d,t), convert(T,x[]), niter))] - x
        return inc
    else
        return zero(x)
    end
end

# short range dispersal kernel, jump to neighbour node
function get_inc(x,D::Nothing,d::DynGraphSpace{T},t) where {T}
    return last(randomwalk(get_graph(d,t),x,2)) - x
end

#increment the trajectory of trait 1 
# such that it follows a reflected brownian motion (1D)
function _reflect1D(x, inc, s::AbstractSegment)
    @assert length(x) == 1 "Only 1D spaces supported"
    _x = x[]
    _inc = inc[] 
    if _x + _inc < s.s
        _inc = 2 * ( s.s - _x ) - _inc
    elseif  _x + _inc > s.e
        _inc = 2 * ( s.e - _x ) - _inc
    else
        return [_inc]
    end
    _reflect1D(x,[_inc],s)
end

function _get_types_dim(s::S) where S<:AbstractSpacesTuple
    T = eltype.(s)
    TT = collect(T)
    _nd = ndims.(s)
    for (i,n) in enumerate(_nd)
        if n > 1
            TT[i] = Vector{TT[i]}
        end
    end
    return TT
end
