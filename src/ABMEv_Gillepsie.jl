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

function update_afterbirth_std!(world,idx_offspring,p::Dict) where T
    # updating competition only the two columns corresponding to agent idx
    α = p["alpha"];K=p["K"];
    x_offspring = get_x(world[idx_offspring])
    for a in skipmissing(world)
        a.d += α(get_x(a),x_offspring)
    end
    # Now updating new agent
    world[idx_offspring].d = sum(α.(get_x.(skipmissing(world)),Ref(x_offspring)))
    world[idx_offspring].b = K(x_offspring)
end

function update_afterdeath_std!(world,x_death,p::Dict) where T
    α = p["alpha"]
    # updating death rate only the two columns corresponding to agent idx
    for a in skipmissing(world)
        a.d -= α(get_x(a),x_death)
    end
end

"""
    function updateWorld_G!(world,C,tspan,p)
Updating rule for gillepsie setting.
Returning dt drawn from an exponential distribution with parameter the total rates of events.
"""
function updateWorld_G!(world,p,update_rates!,tspan,reflected)
    # total update world
    world_alive = skipmissing(world)
    idx_world = collect(eachindex(world_alive))
    # Total rate of events
    U = sum(get_d.(world_alive) .+ get_b.(world_alive))
    dt = - log(rand(Float64))/U
    if dt > 0.
        events_weights = ProbabilityWeights(vcat(get_d.(world_alive),get_b.(world_alive)))
        i_event = sample(events_weights)
        # This is ok since length is a multiple of 2
        N = length(world) - count(ismissing,world)
        if i_event <= N
            # DEATH EVENT
            # In this case i_event is also the index of the individual to die in the world_alive
            idx_offspring = idx_world[i_event]
            x_death = get_x(world[idx_offspring])
            update_afterdeath_std!(world,x_death,p)
            world[idx_offspring] = missing
        else
            # birth event
            idx_offspring = findfirst(ismissing,world)
            # i_event - N is also the index of the individual to give birth in the world_alive
            mum = world[idx_world[i_event-N]]
            world[idx_offspring] = give_birth(mum,p,reflected)
            update_afterbirth_std!(world,idx_offspring,p)
        end
        return dt
    else
        return -1.
    end
end
