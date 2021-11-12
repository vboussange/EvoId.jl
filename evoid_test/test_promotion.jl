#=
Here we show that it is 3 times faster to promote all the traits of individuals to similar types!
=#

using BenchmarkTools

mutable struct mAgent
    # history of traits for geotraits
     x_history
     # birth time of ancestors
     t_history
     # death rate
     d
     #birth rate
     b
end

world_struct = [mAgent(randn(),0.,0.,0.) for i in 1:1000]

function update_world_struct!(world)
    x = rand();
    j = rand(1:length(world))
    if x < 0.5
        # death event
        deleteat!(world,j);
    else
        # birth event
        push!(world, mAgent(world[j].x_history,0.,0.,0.));
    end
end

world_struct = [mAgent(Union{Vector{Float32},Int64}[randn(Float32, 10),rand(1:10)],0.,0.,0.) for i in 1:1000]
@time for i in 1:100000
    update_world_struct!(world_struct)
end
# 0.018369 seconds (77.11 k allocations: 3.963 MiB, 67.32% compilation time)

world_struct = [mAgent(Vector{Float32}[randn(Float32, 10),[rand(1:10)]],0.,0.,0.) for i in 1:1000]
@time for i in 1:100000
    update_world_struct!(world_struct)
end
# 0.006552 seconds (49.95 k allocations: 2.333 MiB)