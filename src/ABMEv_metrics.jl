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
"""
    var(world::Array{Agent};trait=:)
If trait = 0, returns the variance of the geotrait,
knowing that by default it is associated with position trait 1.
If trait > 0, returns the covariance matrix, with first row geotrait and second row trait.
# Arguments

"""
function var(world::Array{Agent{T}};trait=1) where T
    xarray = get_xarray(world,trait)
    return var(xarray,dims=1)
end
"""
    covgeo(world::Array{Agent,1},trait = 0)
If trait = 0, returns the variance of the geotrait,
knowing that by default it is associated with position trait 1.
If trait > 0, returns the covariance matrix, with first row geotrait and second row trait
# Arguments

"""
function covgeo(world::Array{Agent{T},1},trait = 0) where T
    xarray = get_geo.(world)
    if trait > 0
        xstd = reshape(get_x.(world,trait),size(world,1),size(world,2))
        xarray = hcat(xarray,xstd)
    end
    return cov(xarray)
end

"""
    function hamming(world::Array{Agent,1})
Returns a matrix H where H_ij = hamming(a_i,a_j).
The hamming distance is taken through the whole history
and functional space of the agents.
"""
function hamming(world::Array{Agent{T},1}) where T <: Int
    N = size(world,1)
    H = zeros(N,N)
    for (i,a) in enumerate(world)
            for (j,b) in enumerate(world)
                    H[i,j] = hamming(get_xhist(a),get_xhist(b))
            end
    end
    return H
end
