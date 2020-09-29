"""
    abstract type AbstractSpace{Dim,T,F} end
`Dim` is the dimension of the space, `T` is the element type, `ife` is a bool which is `true`
when space is finite
"""

abstract type IsFinite{T} end

#ife stands for is finite
abstract type AbstractSpace{Dim,T,I} end
AbstractSpacesTuple = Tuple{Vararg{AbstractSpace}}
import Base:ndims,isfinite,eltype
Base.ndims(x::AbstractSpace{Dim,T,I}) where {Dim,T,I} = Dim
Base.isfinite(x::AbstractSpace{Dim,T,IsFinite{t}}) where {Dim,T,t} = t #not sure we need this
Base.eltype(::AbstractSpace{Dim,T,I}) where {Dim,T,I} = Dim > 1 ? Tuple{Vararg{T,Dim}} : T
Base.ndims(ss::AbstractSpacesTuple) = length(ss)
Base.eltype(ss::AbstractSpacesTuple) where {Dim,T,I} = Tuple{eltype.(ss)...}

SpaceType=Union{Nothing, AbstractSpace} # not sure what is this used for

"""
$(TYPEDEF)
"""
struct GraphSpace{T} <: AbstractSpace{1,T,IsFinite{true}}
    g::AbstractGraph{T}
end

abstract type AbstractSegment{T<:Number}  <: AbstractSpace{1,T,IsFinite{true}} end

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

A real pace with dimension N and type T
"""
struct RealSpace{N,T} <: AbstractSpace{N,T,IsFinite{false}} end

# TODO:
"""
$(SIGNATURES)

Returns increment
"""
# TODO: find a way to put a type on D in get_inc
function get_inc(D,s::AbstractSpace{Dim,T,I}) where {Dim,T<:AbstractFloat,I<:IsFinite{false}}
    if Dim > 1
        return Tuple(D .* randn(T,Dim))
    else
        return D * randn(T)
    end
end
get_inc(x,D,s::AbstractSpace{Dim,T,I}) where {Dim,T,I<:IsFinite{false}} = get_inc(D,s)

#TODO: there is probably a better way of dealing with those two functions
function get_inc(x,D,s::ContinuousSegment{T}) where {T}
    inc = D * randn(T)
    return _reflect1D(x,inc,s)
end

function get_inc(x,D,s::DiscreteSegment{T}) where {T}
    inc = D * randn()
    return round(T,_reflect1D(x,inc,s))
end

function get_inc(x,D,s::GraphSpace{T}) where {T}
    niter = round(Int,abs(D*randn()))
    if niter > 0
        return last(randomwalk(s.g,x,niter)) - x
    else
        return 0
    end
end

"""
function get_inc_reflected(x::Number,inc::Number,s=-1,e=1)
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
