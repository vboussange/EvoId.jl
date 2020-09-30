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

# This is all about interpolations
import Interpolations:interpolate,Gridded,Linear
struct DiversityFunction
    x
    y
    itp
end
"""
    function geomsmooth2D(xa,ya,itp,smooth)
Return xas,yas,A where A is an interpolated matrix with geometric smooth,
whose axis are xas, yas
# ARGS
`xa` xaxis values, `ya` yaxis values,  `itp` interpolation function, `smooth` smoothing function
"""
function geomsmooth2D(xa,ya,itp,smooth)
    prod(isodd.(smooth)) ? nothing : throw(ArgumentError("smoothing coefficients need to be odd"))
    idx1 = Int((smooth[1]-1)/2)+1:length(xa)- Int((smooth[1]-1)/2)
    idx2 = Int((smooth[2]-1)/2)+1:length(ya)- Int((smooth[2]-1)/2)
    A = [prod(itp.(xa[i-smooth[1] + 1:i]', ya[j-smooth[2] + 1:j])).^(1/prod(smooth)) for i in smooth[1]:length(xa), j in smooth[2]:length(ya)]
    xas = [xa[i] for i in idx1, j in idx2]; yas = [ya[j] for i in idx1,j in idx2]
    return xas,yas,A
end

"""
    function arithsmooth2D(xa,ya,itp,smooth)
Return xas,yas,A where A is an interpolated matrix with arithmetic smooth,
whose axis are xas, yas
# ARGS
`xa` xaxis values, `ya` yaxis values,  `itp` interpolation function, `smooth` smoothing function
"""
function arithsmooth2D(xa,ya,itp,smooth)
    prod(isodd.(smooth)) ? nothing : throw(ArgumentError("smoothing coefficients need to be odd"))
    idx1 = Int((smooth[1]-1)/2)+1:length(xa)- Int((smooth[1]-1)/2)
    idx2 = Int((smooth[2]-1)/2)+1:length(ya)- Int((smooth[2]-1)/2)
    A = [sum(itp.(xa[i-smooth[1] + 1:i]', ya[j-smooth[2] + 1:j]))./prod(smooth) for i in smooth[1]:length(xa), j in smooth[2]:length(ya)]
    xas = [xa[i] for i in idx1,j in idx2]; yas = [ya[j] for i in idx1,j in idx2]
    return xas,yas,A
end

"""
    function interpolate_df(df,xlab,ylab,zlab)
returns an interpolated function itp(x,y) -> z, as well as its axis `xa` and `ya`
"""
function interpolate_df(df,xlab,ylab,zlab)
    sort!(df,[ylab,xlab])
    xa = unique(df[xlab]); ya = unique(df[ylab])
    A = reshape(df[zlab],length(xa),length(ya))
    return DiversityFunction(xa,ya,interpolate((xa,ya),A,Gridded(Linear())))
end

import DataFrames.groupby
"""
    function groupby(f, list::Array)
returns a dictionary that group `list` elements by value of function `f`
"""
groupby(f, list::Array) = begin
  groups = Dict()
  for v in list
    push!(get!(groups, f(v), []), v)
  end
  groups
end
