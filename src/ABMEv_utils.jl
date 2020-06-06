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

"""
    function geomsmooth(x,smooth)
Geometric smoothing, cf `https://en.wikipedia.org/wiki/Exponential_smoothing`

"""
function geomsmooth(x,smooth)
    return [prod(x[i-smooth + 1:i])^(1/smooth) for i in smooth:length(x)]
end

"""
    function arithsmooth(x,smooth)
arithmetic smoothing

"""
function arithsmooth(x,smooth)
    return [sum(x[i-smooth+1:i])/smooth for i in smooth:length(x)]
end

import Plots:cgrad
# asymmetry towards red, blue is only a fifth of the color range
const eth_grad_small = cgrad([colorant"#1F407A", RGB(0.671,0.851,0.914),RGB(1.0,1.0,0.749), RGB(0.992,0.682,0.38),RGB(0.647,0.0,0.149),],[.0,.1])
# symmetry between red and blue
const eth_grad_std = cgrad([colorant"#1F407A", RGB(0.671,0.851,0.914),RGB(1.0,1.0,0.749), RGB(0.992,0.682,0.38),RGB(0.647,0.0,0.149),],[.0,1.])
