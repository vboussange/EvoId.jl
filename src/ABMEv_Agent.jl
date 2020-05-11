abstract type StdAgent end
abstract type MixedAgent end

mutable struct Agent{T,U}
    # history of traits for geotraits
    x_history::Array{U}
    # birth time of ancestors
    t_history::Array{U,1}
    # death rate
    d::Float64
    #birth rate
    b::Float64
end

# Constructors
# This  constructor should be used when one wants to impose the type of the agent (e.g. Mixed)
Agent{T}(xhist::Array{U}) where {T,U} = Agent{T,U}(reshape(xhist,:,1),[0.],0.,1.)

# This constructor is used by default
Agent(xhist::Array{U}) where {U <: Number} = Agent{StdAgent}(xhist)

Agent() = Agent(Float64[],0.,0.,1.)
import Base.copy
copy(a::Agent{T,U}) where {T,U} = Agent{T,U}(copy(a.x_history),copy(a.t_history),copy(a.d),copy(a.b))
copy(m::Missing) = missing

"""
    function new_world_G(nagents::Int,p::Dict; spread = 1., offset = 0.)
Returns an array of type Array{Union{Missing,Agent}} initialised with normal distribution.
Only relevant for Gillepsie algorithm as of now.
"""
function new_world_G(nagents::Int,p::Dict; spread = 1., offset = 0.)
    typeof(spread) <: Array ? spread = spread[:] : nothing;
    typeof(offset) <: Array ? offset = offset[:] : nothing;
    agent0 = [Agent( spread  .* randn(length(spread)) .+ offset) for i in 1:nagents]
    world0 = vcat(agent0[:],repeat([missing],Int(p["NMax"] - nagents)))
    return world0
end

# returns trait i of the agent
get_x(a::Agent) = a.x_history[:,end]
function get_geo(a::Agent{U,T},t::Number) where {U,T}
    tarray = vcat(a.t_history[2:end],convert(T,t))
    tarray .-= a.t_history
    return sum(get_xhist(a,1) .* tarray)
end
# This method can acces geotrait, while the second not
get_x(a::Agent,t::Number,i::Integer) = i > 0 ? a.x_history[Int(i),end] : get_geo(a,t)
get_x(a::Agent,i::Integer) = a.x_history[Int(i),end]
get_xhist(a::Agent,i::Number) = a.x_history[Int(i),:]
get_xhist(a::Agent) = a.x_history
get_d(a::Agent) = a.d
get_b(a::Agent) = a.b
get_fitness(a::Agent) = a.b - a.d
get_dim(a::Agent) = size(a.x_history,1)
get_nancestors(a::Agent) = size(a.x_history,2)

get_x(world::Array{T},trait::Integer) where {T <: Agent} = trait > 0 ? reshape(hcat(get_x.(world,trait)),size(world,1),size(world,2)) : throw(ErrorException("Not the right method, need `t` as an argument"))
"""
    get_x(world::Array{T},t::Number,trait::Integer) where {T <: Agent}
Returns trait of every agents of world in the form of an array which dimensions corresponds to the input.
If trait = 0 , we return the geotrait.

"""
get_x(world::Array{T},t::Number,trait::Integer) where {T <: Agent} = trait > 0 ? reshape(hcat(get_x.(world,trait)),size(world,1),size(world,2)) : reshape(hcat(get_geo.(world,t)),size(world,1),size(world,2))

"""
    get_xarray(world::Array{Agent,1})
Returns every traits of every agents of world in the form of an array
"""
function get_xarray(world::Array{T,1}) where {T <: Agent}
    return hcat(get_x.(world)...)
end

function get_xarray(world::Array{T,1},t::Number,geotrait::Bool=false) where {T <: Agent}
    xarray = hcat(get_x.(world)...)
    if geotrait
        xarray = vcat( xarray, get_geo.(world,t)')
    end
    return xarray
end

# """
#     get_xhist(world::Vector{Agent},geotrait = false)
# Returns the trait history of every agents of world in the form of an 3 dimensional array,
# with
# - first dimension as the agent index
# - second as time index
# - third as trait index
# If geotrait = true, then a last trait dimension is added, corresponding to geotrait.
# Note that because number of ancestors are different between agents, we return an array which size corresponds to the minimum of agents ancestors,
# and return the last generations, dropping the youngest ones
# """
# function get_xhist(world::Vector{T}) where {T <: Agent}
#     hist = minimum(get_nancestors.(world))
#     ntraits = get_dim(first(world));
#     xhist = zeros(length(world), hist, ntraits + geotrait);
#     for (i,a) in enumerate(world)
#         xhist[i,:,1:end-geotrait] = get_xhist(a)[:,end-hist+1:end]';
#     end
#     return xhist
# end

# TODO: This method broken, when one ask for the geotraits
# function get_xhist(world::Vector{T},t::Number,geotrait = false) where {T <: Agent}
#     hist = minimum(get_nancestors.(world))
#     ntraits = get_dim(first(world));
#     xhist = zeros(length(world), hist, ntraits + geotrait);
#     for (i,a) in enumerate(world)
#         xhist[i,:,1:end-geotrait] = get_xhist(a)[:,end-hist+1:end]';
#         if geotrait
#             xhist[i,:,ntraits+geotrait] = cumsum(get_xhist(a,1))[end-hist+1:end]
#         end
#     end
#     return xhist
# end


function world2df(world::Array{T,1}) where {T <: Agent}
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



"""
    function increment_x!(a::Agent{StdAgent,U},t::U,p::Dict) where U
This function increments agent by random numbers specified in p
ONLY FOR CONTINUOUS DOMAINS
"""
function increment_x!(a::Agent{StdAgent,U},t,p::Dict) where U
    tdim = length(p["D"])
    reflected = haskey(p,"reflected") ? p["reflected"] : false
    if reflected
        inc = [get_inc_reflected(get_x(a,1),p["D"][1] *randn())]
        if  tdim > 1
            inc = vcat(inc,(rand(tdim-1) < p["mu"][2:end]) .* p["D"][2:end] .* randn(tdim-1))
        end
    else
        # inc = yes no mutation * mutation
        inc = (rand(tdim) < vec(p["mu"])) .* vec(p["D"][:]) .* randn(tdim)
    end
    a.x_history = hcat(a.x_history, get_x(a) + reshape(inc,:,1));
    push!(a.t_history,t)
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
get_inc_reflected(x::Number,inc::Number,s=-1,e=1)
    Here we increment the trajectory of trait 1 such that it follows a reflected brownian motion (1D)
"""
function get_inc_reflected(x::Number,inc::Number,s=-1,e=1)
    if x + inc < s
        inc = 2 * ( s - x ) - inc
    elseif  x + inc > e
        inc = 2 * ( e - x ) - inc
    else
        return inc
    end
    get_inc_reflected(x,inc,s,e)
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
