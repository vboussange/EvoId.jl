p = Dict(1:1:1000 .=> 1:1:1000)
delete!(p,999)
a = Any[rand() for i in 1:1000]
a[2:2:1000] .= missing
get_idx(w) = collect(eachindex(skipmissing(w)))
using BenchmarkTools
@btime a[get_idx(a)[3]]
@btime p[3]


@btime [true for i in 1:1000]
[true for i in 1:1000]

@btime sort!(values(p))
function findfreeidx(p)
    i = 1
    while haskey(p,i)
        i+=1
    end
    i
end

function get_xp(a,j)
    for (i,a) in enumerate(skipmissing(a))
        if i==j
            return a
            break
        end
    end
end
@btime get_xp(a,998)

@btime findfreeidx(p)
@btime collect(keys(p))[2]

@btime(count(ismissing,))
