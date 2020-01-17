# In this file lie the function for Gillepsie algorithm


"""
    give_birth(a::Agent,p)
Used for Gillepsie setting
"""
function give_birth(a::Agent,p::Dict,reflected)
    new_a = copy(a)
    increment_x!(new_a,p,reflected=reflected)
    return new_a
end

function update_afterbirth_std!(world,C,idx::Int,p::Dict) where T
    traits = get_x.(world)
    # updating competition only the two columns corresponding to agent idx
    for i in 1:length(traits)
        C[i,idx] = α(traits[i],traits[idx],p["n_alpha"],p["sigma_a"])
        C[idx,i] = C[i,idx]
    end
    # updating death rate only the two columns corresponding to agent idx
    for (i,a) in enumerate(world)
        a.d += C[i,idx] / p["K0"]
    end
    # Now updating new agent
    world[idx].d = sum(C[idx,:]) / p["K0"]
    world[idx].b = K(traits[idx][:,end],1.,p["n_K"],p["sigma_K"])
end

function update_afterdeath_std!(world,C,idx::Int,p::Dict) where T
    traits = get_x.(world)
    # updating death rate only the two columns corresponding to agent idx
    for (i,a) in enumerate(world)
        a.d -= C[i,idx] / p["K0"]
    end
    # updating competition only the two columns corresponding to agent idx
    for i in 1:length(traits)
        C[i,idx] = .0
        C[idx,i] = .0
    end
end

"""
    function updateWorld_G!(world,C,tspan,p)
Updating rule for gillepsie setting.
Returning dt drawn from an exponential distribution with parameter the total rates of events.
"""
function updateWorld_G!(world,C,p,update_rates!,tspan,reflected)
    # total update world
    world_alive = skipmissing(world)
    idx_world = collect(eachindex(world_alive))
    # Total rate of events
    U = sum(get_d.(world_alive) .+ get_b.(world_alive))
    dt = - log(rand(Float64))/U
    events_weights = ProbabilityWeights(vcat(get_d.(world_alive),get_b.(world_alive)))
    i_event = sample(events_weights)
    # This is ok since length is a multiple of 2
    I = Int(length(events_weights) / 2)
    if i_event <= I
        # DEATH EVENT
        idx_offspring = idx_world[i_event]
        world[idx_offspring] = missing
        update_afterdeath_std!(world_alive,C,idx_offspring,p)
    else
        # birth event
        idx_offspring = findfirst(ismissing,world)
        world[idx_offspring] = give_birth(world[idx_world[i_event-I]],p,reflected)
        update_afterbirth_std!(world_alive,C,idx_offspring,p)

    end
    return dt
end
