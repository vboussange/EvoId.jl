"""
        generalised_gaussian(x::Float64,mu::Float64,sigma::Float64,epsilon::Float64)
"""
function generalised_gaussian(x::Float64,mu::Float64,sigma::Float64,epsilon::Float64)
        return exp( -.5 * ((x-mu) / sigma)^epsilon)
end

"""
        gaussian(x::Float64,mu::Float64,sigma::Float64) = generalised_gaussian(x,mu,sigma,2)
"""
gaussian(x::Float64,mu::Float64,sigma::Float64) = generalised_gaussian(x,mu,sigma,2.)
