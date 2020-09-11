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

import Statistics.var
# TODO: rename this to gamma diversity
"""
    var(world::Array{Agent};trait=:)
If trait = 0, returns the variance of the geotrait,
knowing that by default it is associated with position trait 1.
If trait > 0, returns the covariance matrix, with first row geotrait and second row trait.
# Arguments

"""
function var(world::Array{U,1};trait=1) where U <: Union{Missing,Agent{T}} where T
    world = collect(skipmissing(world))
    xarray = get_x(world,trait)
    return var(xarray,dims=1,corrected=false)
end
"""
    covgeo(world::Array{Agent,1},trait = 0)
If trait = 0, returns the variance of the geotrait,
knowing that by default it is associated with position trait 1.
If trait > 0, returns the covariance matrix, with first row geotrait and second row trait
# Arguments

"""
function covgeo(world::Array{U,1},t::Number,trait = 0) where U <: Union{Missing,Agent}
    world = collect(skipmissing(world))
    xarray = Float64.(get_geo.(world,t))
    if trait > 0
        xstd = reshape(Float64.(get_x.(world,trait)),size(world,1),size(world,2))
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
function hamming(world::Array{Agent,1}) where T <: Int
    N = size(world,1)
    H = zeros(N,N)
    for (i,a) in enumerate(world)
            for (j,b) in enumerate(world)
                    H[i,j] = hamming(get_xhist(a),get_xhist(b))
            end
    end
    return H
end
"""
    get_alpha_div(world::Array{U,1},t::Number,trait=1) where U <: Union{Missing,Agent}
Mean of the local variance of `trait` per patch
# Arguments
"""
function get_alpha_div(world::Array{U,1},t::Number,trait=1) where U <: Union{Missing,Agent}
    _xall_df = world2df(world,t,true)
    xall_per_patch = groupby(_xall_df, :x1,sort=true)
    if trait == 0
        # need to convert to Float64, otherwise infinite variance
        return mean([var(Float64.(xp.g),corrected=false) for xp in xall_per_patch])
    else
        return mean([var(Float64.(xp[:,trait+1]),corrected=false) for xp in xall_per_patch])
    end
end

"""
    get_beta_div(world::Array{U,1},t::Number,trait=1) where U <: Union{Missing,Agent}
Variance of the mean of `trait` per patch
# Arguments
"""
function get_beta_div(world::Array{U,1},t::Number,trait=1) where U <: Union{Missing,Agent}
    _xall_df = world2df(world,t,true)
    xall_per_patch = groupby(_xall_df, :x1,sort=true)
    if trait == 0
        # need to convert to Float64, otherwise infinite variance
        sbar_i = [mean(Float64.(xp.g)) for xp in xall_per_patch]
    else
        sbar_i = [mean(Float64.(xp[:,trait+1])) for xp in xall_per_patch]
    end
    return var(sbar_i,corrected=false)
end

"""
    function get_hamming_dist_hist(a1,a2,time,trait=1)
Returns a metric corresponding to the exact time agents have been isolated to each other
"""
function get_hamming_dist_hist(a1,a2,trait=1,time = 0)
        thist = vcat(get_thist(a1),get_thist(a2))
        ttot = sort!(unique(thist))
        xhist = zeros(2,length(ttot))
        for (i,a) in enumerate([a1,a2])
                _thist = get_thist(a)
                _xhist = get_xhist(a,trait)
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
        return sum((xhist[1,1:end-tt] .!= xhist[2,1:end-tt])[:] .* (ttot[2:end] .- ttot[1:end-1]))
        # return xhist
        # return (thist[2:end] .- thist[1:end-1])
end

"""
    function get_pairwise_average_isolation(world,trait=1,time = 0)
Returns the average pairwise isolation by time between all agents of `world`,
using metrics `get_hamming_dist_hist`
"""
function get_pairwise_average_isolation(world,trait=1,time = 0)
        N = size(world,1)
        D = zeros(N)
        # Here you should do a shared array to compute in parallel
        Threads.@threads for i in 1:(N-1)
                for j in i+1:N
                            C = get_hamming_dist_hist(world[i],world[j],trait,time)
                            D[i] += C
                            D[j] += C
                end
        end
        return sum(D) / N^2
end

"""
    function get_local_pairwise_average_isolation(world,trait=1,time = 0)
Similar to `get_pairwise_average_isolation`, but the pairwise distance is calculated within demes.
An average of this metrics by deme is return.
"""
function get_local_pairwise_average_isolation(world,trait=1,time = 0)
        f(a) = get_x(a,1)
        groups = groupby(f,world)
        d = get_pairwise_average_isolation.(collect(values(groups)),trait,time)
        mean(d)
end
