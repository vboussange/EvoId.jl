#=
Here we test that having a vector of agents that changes through time by push! or deletat! 
is better than a fixed vector with Union{nothing, Agent}
The test is positive
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

world_struct_nothing = Union{mAgent,Nothing}[[mAgent(randn(),0.,0.,0.) for i in 1:1000];fill(nothing,1000)]
world_struct = [mAgent(randn(),0.,0.,0.) for i in 1:1000]

function update_world_struct_nothing!(world)
    x = rand();
    j = rand((1:length(world))[.!isnothing.(world)])
    if x < 0.5
        # death event
        world[j] = nothing;
    else
        # birth event
        jj = first((1:length(world))[isnothing.(world)])
        world[jj] = mAgent(randn(),0.,0.,0.);
    end
end
update_world_struct_nothing!(world_struct_nothing)

function update_world_struct!(world)
    x = rand();
    j = rand(1:length(world))
    if x < 0.5
        # death event
        deleteat!(world,j);
    else
        # birth event
        push!(world,  mAgent(randn(),0.,0.,0.));
    end
end
update_world_struct!(world_struct)

world_struct_nothing = Union{mAgent,Nothing}[[mAgent(randn(),0.,0.,0.) for i in 1:1000];fill(nothing,1000)]
@time for i in 1:100000
    update_world_struct_nothing!(world_struct_nothing)
end
# 1.046872 seconds (700.19 k allocations: 1.811 GiB, 27.58% gc time)

world_struct = [mAgent(randn(),0.,0.,0.) for i in 1:1000]
@time for i in 1:100000
    update_world_struct!(world_struct)
end
# 0.008058 seconds (100.29 k allocations: 3.106 MiB)