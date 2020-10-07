import StatsBase:Histogram,fit,step,normalize

d1(x,y) = norm(x .- y,1)
d2(x,y) = norm(x .- y,2)
hamming(k::Array, l::Array) = sum(k .!= l)

"""
    H_discrete(s)
Interconnectedness measure as in Nordbotten 2018 for discrete setup
"""
function H_discrete(s)
    N = length(s)
    H = 0
    for x in s
        for y in s
            H += d2(x,y)
        end
    end
    return H / N^2
end


"""
    findclusters(v::Vector,allextrema =true)
Returns a tuple with the cluster mean and its associated weight
# Arguments

"""
function findclusters(v::Vector,allextrema =true)
    # this function fits a histogram to some distribution and extracts its local maxima
    h = fit(Histogram,v)
    s = step(h.edges...)
    x = normalize(h.weights,1)
    dx = x[2:end] .- x[1:end-1]
    if allextrema
        sdx = dx[2:end] .* dx[1:end-1]
        idx = findall(i -> i < 0, sdx)
    else
        idx = findall(i -> dx[i] > 0 && dx[i+1] < 0, 1:length(dx))
    end
    return collect(h.edges...)[idx .+ 1] .+ s/2, x[idx .+ 1]
end

import Statistics:var,mean
# TODO: rename this to gamma diversity
"""
    var(world::Array{Agent};trait=:)
If trait = 0, returns the variance of the geotrait,
knowing that by default it is associated with position trait 1.
If trait > 0, returns the covariance matrix, with first row geotrait and second row trait.
# Arguments

"""
function var(world::World;trait=1)
    xarray = Float64.(get_x(world,trait))
    return var(xarray,dims=1,corrected=false)
end

"""
function mean(world::World;trait=1)

"""
function mean(world::World;trait=1)
    xarray = Float64.(get_x(world,trait))
    return mean(xarray,dims=1)
end
"""
    covgeo(world::Array{Agent,1},trait = 0)
If trait = 0, returns the variance of the geotrait,
knowing that by default it is associated with position trait 1.
If trait > 0, returns the covariance matrix, with first row geotrait and second row trait
# Arguments

"""
function covgeo(world::World,trait = 0)
    xarray = Float64.(get_geo(world))
    if trait > 0
        xstd = get_x(world,trait)
        xarray = hcat(xarray,xstd)
    end
    return cov(xarray,corrected=false)
end

"""
    function hamming(world::Array{Agent,1})
Returns a matrix H where H_ij = hamming(a_i,a_j).
The hamming distance is taken through the whole history
and functional space of the agents.
"""
function hamming(world::World) where T <: Int
    n = size(world)
    H = zeros(n,n)
    for (i,a) in enumerate(agents(world))
            for (j,b) in enumerate(world)
                    H[i,j] = hamming(get_xhist(a),get_xhist(b))
            end
    end
    return H
end
"""
    get_alpha_div(world::Array{U,1},t::Number,trait=1) where U <: Union{Missing,Agent}
Mean of the local variance of `trait` per patch. If trait=0, we get the geotrait
# Arguments
"""
function get_alpha_div(world::World,trait=1)
    g = groupby(a->a[1],agents(world))
    if trait == 0
        return mean([var(Float64.(get_geo(World(subw,space(world),parameters(world)))),corrected=false) for subw in values(g)])
    else
        # here the second mean is here when subspace is multidimensional
        v = [var(World(subw,space(world),parameters(world)),trait=trait) for subw in values(g)]
        h = vcat(v...)
        return mean(h)
    end
end

"""
    get_beta_div(world::Array{U,1},t::Number,trait=1) where U <: Union{Missing,Agent}
Variance of the mean of `trait` per patch
# Arguments
"""
function get_beta_div(world::World,trait=1)
    g = groupby(a->a[1],agents(world))
    if trait == 0
        # need to convert to Float64, otherwise infinite variance
        sbar_i = [mean(Float64.(get_geo(World(subw,space(world),parameters(world))))) for subw in values(g)]
        return var(sbar_i,corrected=false)
    else
        m = [mean(World(subw,space(world),parameters(world)),trait=trait) for subw in values(g)]
        h=vcat(m...)
        return mean(var(h,dims=1,corrected=false))
    end
end

"""
$(SIGNATURES)
returns `xhist,ttot`, where `xhist` is a matrix with dimension `lenght(world)` x `length(thist)+1`,
which consists of geographical history of ancestors at every time step.
If `time` is specified and is higher that the highest time step in world,
then the last column of xhist corresponds to actual position of agents
"""

function get_xhist_mat(agentarray::Vector{A},trait=1,time = 0) where {A<:AbstractAgent}
        thist = vcat(get_thist.(agentarray)...)
        ttot = sort!(unique(thist))
        xhist = zeros(length(agentarray),length(ttot))
        for (i,a) in enumerate(agentarray)
            _thist = get_thist(a)
            # Here we check that there is no redundancy of thist because of casting errors
            # In the future, we should remove this check, as the type of time has been set to Float64
            # Hence there should be no more problems of this type
            ttrue = _thist[2:end] .- _thist[1:end-1] .> 0
            if count(ttrue) < length(_thist) - 1
                _tt = vcat(ttrue,true)
                _thist = _thist[_tt] # we drop the position that was occupied for dt = .0
                _xhist = get_xhist(a,trait)[_tt]
            else
                _xhist = get_xhist(a,trait)
            end
            k = 0
            _l = length(_thist)
            for j in 1:(_l - 1)
                    xhist[i,j+k] = _xhist[j]
                            while _thist[j+1] > ttot[j+k+1]
                                    k+=1
                                    xhist[i,j+k] = _xhist[j]
                                    if j+k == length(ttot)
                                            break
                                    end
                            end
            end
            xhist[i,_l+k:end] .= _xhist[end]
        end
        tt = 0
        time > ttot[end] ? ttot = vcat(ttot,time) : tt = 1
        return xhist[:,1:end-tt],ttot
end

"""
    function get_dist_hist(a1,a2,dist,trait=1,time = 0)

Returns the integral of the distance `dist` through time of `trait` between `a1.x` and `a2.x`.

```math
\\int d(x_1,x_2,t)dt
```
"""
function get_dist_hist(a1,a2,dist,trait=1,time = 0)
        xhist,ttot = get_xhist_mat([a1,a2],trait,time)
        return sum(dist.(xhist[1,:],xhist[2,:]) .* (ttot[2:end] .- ttot[1:end-1]))
end

"""
    function truncvar(X::AbstractArray)
Returns the truncated variance of `X`.
"""

function truncvar(X::AbstractArray)
        N = length(X)
        c = [count(y->y!=x,X) for x  in unique(X)]
        return sum(c .* (N .-c)) /(2*N^2)
end

"""
    function get_pairwise_average_isolation(world;trait=1,trunc=false)
Returns the integrated pairwise squared distance between all agents of `world` wrt `trait`.
If `trunc=true` then the distance is truncated to binary value 0 or 1.
"""
function get_pairwise_average_isolation(world::World;trait=1,trunc=false)
        xhist,ttot = get_xhist_mat(agents(world))
        if trunc
                V = truncvar.([xhist[:,i] for i in 1:size(xhist,2)])
        else
                V = var(xhist,dims=1,corrected=false)[:]
        end
        return sum(V .* (ttot[2:end] .- ttot[1:end-1]))
end

"""
    function get_local_pairwise_average_isolation(world,dist,trait=1)
Similar to `get_pairwise_average_isolation`, but the pairwise distance is calculated within demes.
An average of this metrics by deme is return.
"""
function get_local_pairwise_average_isolation(world::World;trait=1,trunc=false)
        f(a) = a[1]
        groups = groupby(f,agents(world))
        smallworlds = []
        for v in values(groups)
            myagents = [a for a in v]
            push!(smallworlds,World(myagents,world.space,world.p,world.t))
        end
        d = get_pairwise_average_isolation.(smallworlds,trait=trait,trunc=trunc)
        mean(d)
end
