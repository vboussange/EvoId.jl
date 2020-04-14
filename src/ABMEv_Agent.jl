mutable struct Agent{T}
    # history of traits for geotraits
    x_history::Array{T}
    # death rate
    d::Float64
    #birth rate
    b::Float64
end

# Constructors

Agent(xhist::Array{T}) where T <: Number = Agent{T}(reshape(xhist,:,1),0.,1.)
Agent() = Agent(Float64[],0.,1.)
import Base.copy
copy(a::Agent) = Agent(a.x_history,a.d,a.b)
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
get_x(a::Agent,i::Number) = a.x_history[Int(i),end]
get_x(a::Agent) = a.x_history[:,end]
get_xhist(a::Agent,i::Number) = a.x_history[Int(i),:]
get_xhist(a::Agent) = a.x_history
get_geo(a::Agent) = sum(get_xhist(a,1))
get_d(a::Agent) = a.d
get_b(a::Agent) = a.b
get_fitness(a::Agent) = a.b - a.d
get_dim(a::Agent) = size(a.x_history,1)
get_nancestors(a::Agent) = size(a.x_history,2)
"""
    get_xarray(world::Array{Agent{T}},trait::Int) where T
Returns trait of every agents of world in the form of an array which dimensions corresponds to the input.
Particularly suited for an array world corresponding to a timeseries.

"""
get_xarray(world::Array{Agent{T}},trait::Int) where T = reshape(hcat(get_x.(world,trait)),size(world,1),size(world,2))

"""
    get_xhist(world::Vector{Agent{T}},geotrait = false)
Returns the trait history of every agents of world in the form of an 3 dimensional array,
with
- first dimension as the agent index
- second as time index
- third as trait index
If geotrait = true, then a last trait dimension is added, corresponding to geotrait.
Note that because number of ancestors are different between agents, we return an array which size corresponds to the minimum of agents ancestors,
and return the last generations, dropping the youngest ones
"""
function get_xhist(world::Vector{Agent{T}},geotrait = false) where T
    hist = minimum(get_nancestors.(world))
    ntraits = get_dim(first(world));
    xhist = zeros(length(world), hist, ntraits + geotrait);
    for (i,a) in enumerate(world)
        xhist[i,:,1:end-geotrait] = get_xhist(a)[:,end-hist+1:end]';
        if geotrait
            xhist[i,:,ntraits+geotrait] = cumsum(get_xhist(a,1))[end-hist+1:end]
        end
    end
    return xhist
end


"""
increment_x!(a::Agent{Float64},p::Dict)
    This function increments current position by inc and updates xhist,
    ONLY FOR CONTINUOUS DOMAINS
"""
function increment_x!(a::Agent{Float64},p::Dict)
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
 end

 """
     function increment_x!(a::Agent{Int64},inc::Array{Float64})
 This function increments current position by inc and updates xhist,
     ONLY FOR GRAPH TYPE DOMAINS
 """
 function increment_x!(a::Agent{Int64},p::Dict)
     # we have to add 1 otherwise it stays on the same edge
     inc = randomwalk(p["g"],get_x(a,1),Int(p["D"][1]) + 1)
     # the last coef of inc corresponds to last edge reached
     a.x_history = hcat(a.x_history, reshape(inc[end:end],:,1));
  end


"""
get_inc_reflected(x::Float64,inc::Float64,s=-1,e=1)
    Here we increment the trajectory of trait 1 such that it follows a reflected brownian motion (1D)
"""
function get_inc_reflected(x::Float64,inc::Float64,s=-1,e=1)
    if x + inc < s
        inc = 2 * ( s - x ) - inc
    elseif  x + inc > e
        inc = 2 * ( e - x ) - inc
    else
        return inc
    end
    get_inc_reflected(x,inc,s,e)
end

# need to make sure that this is working correctly
"""
    α(a1::Array{Float64},a2::Array{Float64},n_alpha::Float64,sigma_a::Array{Float64})
Generalised gaussian competition kernel
"""
function α(a1::Array{Float64},a2::Array{Float64},n_alpha::Float64,sigma_a::Array{Float64})
        return exp( -.5* sum(sum((a1 .- a2).^n_alpha,dims=2)./ sigma_a[:].^n_alpha))
end
"""
    α(a1::Array{Float64},a2::Array{Float64},n_alpha::Float64,sigma_a::Array{Float64})
Default Gaussian competition kernel
"""
α(a1::Array{Float64},a2::Array{Float64},sigma_a::Array{Float64}) = α(a1,a2,2,sigma_a)

"""
    K(x::Array{Float64},K0::Float64,μ::Array{Float64},sigma_K::Array{Float64})
Gaussian resource kernel
"""
function K(x::Array{Float64},K0::Float64,μ::Array{Float64},sigma_K::Array{Float64})
    # return K0*exp(-sum(sum((x .- μ).^n_K,dims=2)./sigma_K[:].^n_K))
    return K0 * pdf(MvNormal(μ,sigma_K),x)
end
"""
    K(x::Array{Float64},K0::Float64,sigma_K::Array{Float64})
Gaussian resource kernel with mean 0.
"""
function K(x::Array{Float64},K0::Float64,sigma_K::Array{Float64})
    # return K0*exp(-sum(sum((x .- μ).^n_K,dims=2)./sigma_K[:].^n_K))
    return K0 * pdf(MvNormal(sigma_K),x)
end

KK(x::Array{Float64},K0::Float64,n_K::Float64,sigma_K::Array{Float64},μ1::Float64,μ2::Float64) = K(x,K0/2,μ1,sigma_K) + K(x,K0/2,μ2,sigma_K)

"""
    function tin(t::Float64,a::Float64,b::Float64)
if t in [a,b) returns 1. else returns 0
"""

function tin(t::Float64,a::Float64,b::Float64)
    return t>=a && t<b ? 1. : 0.
end

function split_move(t)
    return .0 + 1/100*(t-20.)*tin(t,20.,120.) + tin(t,120.,Inf64)
end

function split_merge_move(t)
    return .0 + 1/30*(t-10.)*tin(t,10.,40.) + tin(t,40.,70.) + (1- 1/30*(t-70.))*tin(t,70.,100.)
end
