mutable struct World{A<:AbstractAgent, S<:AbstractSpacesTuple,T} 
    agents::Vector{A}
    space::S
    D #increment range
    mu #increment probability
    NMax #maximum number of individuals
    t::T #time
end

"""
World(agents, s, D, mu, NMax; t=0.)

Constructs a world.

# Arguments

* `w` a vector of agents, 
* `s` a tuple of evolutionary spaces, 
* `mu` a vector that contains the increment rate. If `mu[i]`
    has dimension greater than 1, 
    then mutations happen independently at each dimension
    of `s[i]`.
* `D`, a vector that contains the increment ranges. Only `nothing` is 
    supported for `GraphSpace`, equivalent to a random walk of length 1.
* `NMax` the maximum number of individuals allowed during the simulation
    
# Examples

```julia

nodes = 7
g = star_graph(nodes)
landscape = GraphSpace(g)
θ = [rand([-1,1]) for i in 1:nodes]
traitspace = RealSpace(1)
evolspace = (landscape,traitspace)

D = [nothing,5e-2]
mu = [1f-1,1f-1]
p = Dict("NMax" => 2000,
    "D" => D,
    "mu" => mu)
myagents = [Agent(evolspace,[rand(1:nodes),randn() * D[2]]) for i in 1:K]

w0 = World(myagents,evolspace,p)
```

"""
function World(w::Vector{<:AbstractAgent{X}}, s::S, D, mu, NMax; t=0.) where {X, S<:AbstractSpacesTuple}
    # if typeof(p["D"]) != eltype(skipmissing(w)[1])
    #     throw(ArgumentError("Diffusion coefficient does not match with underlying space\n `D::Tuple`"))
    # end

    for _m in mu
        if !(eltype(_m) <: AbstractFloat)
            throw(ArgumentError("elements of `mu` should be of type AbstractFloat\n
                                to decide if mutations occur from a uniform probability law"))
        end
    end

    @assert length(mu) == length(s) "Length of parameter `mu` should correspond to dimension of underlying space"
    @assert length(D) == length(s) "Length of parameter `D` should correspond to dimension of underlying space"
    SS = eltype.(s)
    _SS = _get_types_dim(s)
    # checking that eltypes of D are nothing or abstractfloat
    for (i,_S) in enumerate(SS)
        if _S <: Integer
            # handling for discrete space is different
            Di = D[i]
            @assert (isnothing(Di) || typeof(Di) <: AbstractFloat) "type`D` at dimension $i should be whether nothing or a float"
            _SS[i] = typeof(Di)
        end
    end
    # converting internally D ranges to align with subspaces eltype
    D2 = Union{_SS...}[]
    for (i,Di) in enumerate(D)
        # making sure that we only modify D for continuous space
        # for discrete space, D should be whether Nothing or AbstractFloat
        if SS[i] <: AbstractFloat
            try
                push!(D2, convert(_SS[i],Di))
            catch e
                throw(ArgumentError("`D` at dimension $i does not match with underlying space"))
            end
        else
            push!(D2, Di)
        end
    end
    # TODO: how not having to use this horrible thing?
    World(w, s, D2, mu, NMax, convert(X,t))
end

# this throws an iterators of agents in the world
agents(world::World) = world.agents
parameters(world::World) = world.p
time(w::World) = w.t
space(w::World) = w.space
maxsize(w::World) = w.NMax
get_D(w::World) = w.D
get_mu(w::World) = w.mu

# this throws indices that are occupied by agents
# this throws agents of an abstract array of size size(world)
import Base:size,getindex
Base.length(world::World) = length(world.agents)
Base.copy(w::W) where {W<:World} = W(copy.(w.agents), w.space, copy(w.D), copy(w.mu) ,copy(w.t))
## Accessors
"""
    Base.getindex(w::World,i) 

Get x of world without geotrait.
"""
Base.getindex(w::World,i) = w.agents[i]

function Base.show(io::IO, w::World{A,S,T}) where {A,S,T}
     println(io, "World with agents of type", A)
 end

addAgent!(w::World,a::AbstractAgent) = begin
    push!(w.agents,a)
end
removeAgent!(w::World,i::Int) = begin
    deleteat!(w.agents,i)
end

update_clock!(w::World{A,S,T},dt) where {A,S,T} = begin
    w.t = w.t + dt
    return nothing
end


"""
    get_geo(w)
"""
get_geo(w::World) = map(a-> get_geo(a,time(w)), agents(w))

"""
$(SIGNATURES)

Returns trait of every agents of world in the form of an array which dimensions corresponds to the input.
If `trait = 0` , we return the geotrait.
!!! warning "Warning"
    Geotrait might be deprecated in the future.
"""
function get_x(w::World,trait)
    if !(trait == 0)
        if ndims(space(w)[trait]) > 1
            return hcat(collect.(getindex.(agents(w),trait))...)'
        else
            return collect(getindex.(agents(w),trait))
        end
    else
        return collect(get_geo(w))
    end
end

"""
$(SIGNATURES)

Returns every traits of every agents of `world` in the form **of a one dimensional array** (in contrast to `get_x`).
If `geotrait=true` the geotrait is also added to the set of trait, in the last column.
If you do not want to specify `t` (only useful for geotrait), it is also possible to use `get_xarray(world::Array{T,1}) where {T <: Agent}`.
!!! warning "Warning"
    It does not work with subspace where ndims(subspace) > 1.
"""
function get_xarray(world::World,geotrait::Bool=false)
    xarray = get_x(world,Colon())
    if geotrait
        xarray = hcat(xarray, get_geo.(agents(world),world.t))
    end
    return xarray
end
@deprecate get_xarray(world,geotrait=false) get_x(world,Colon())

"""
    give_birth(mum_idx::Int,w::World)
Copies agent within index `mum_idx`, and increment it by dx.
Return new agent (offspring).
"""
function give_birth(mum_idx::Int,w::World)
    new_a = copyxt(w[mum_idx])
    increment_x!(new_a,space(w),get_D(w), get_mu(w),time(w))
    return new_a
end
