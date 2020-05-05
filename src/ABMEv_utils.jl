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

import DSP.conv
"""
    ma(x::Array{T},f) where T <: Number
Moving average over array x, using f as the filter, i.e. the number of points to average on. Better choosing an odd number
"""
function ma(x::Array{T},f) where T <: Number
    _N = length(x)
    _s = Int((f-1)/2)
    return conv(x,ones(f)./f)[_s:_s+_N-1]
end

import Plots:cgrad
const eth_grad_small = cgrad([colorant"#1F407A", RGB(0.671,0.851,0.914),RGB(1.0,1.0,0.749), RGB(0.992,0.682,0.38),RGB(0.647,0.0,0.149),],[.0,.1])
const eth_grad_std = cgrad([colorant"#1F407A", RGB(0.671,0.851,0.914),RGB(1.0,1.0,0.749), RGB(0.992,0.682,0.38),RGB(0.647,0.0,0.149),],[.0,1.])
