"""
        generalised_gaussian(x::Number,mu::Number,sigma::Number,epsilon::Number)
"""
function generalised_gaussian(x::Number,mu::Number,sigma::Number,epsilon::Number)
        return exp( -.5 * ((x-mu) / sigma)^epsilon)
end

"""
        gaussian(x::Number,mu::Number,sigma::Number) = generalised_gaussian(x,mu,sigma,2)
"""
gaussian(x::Number,mu::Number,sigma::Number) = generalised_gaussian(x,mu,sigma,2.)

"""
    ma(x::Array{T},f) where T <: Number
Moving average over array x, using f as the filter, i.e. the number of points to average on. Better choosing an odd number
"""
function ma(x::Array{T},f) where T <: Number
    _N = length(x)
    _s = Int((f-1)/2)
    return conv(x,ones(f)./f)[_s:_s+_N-1]
end
