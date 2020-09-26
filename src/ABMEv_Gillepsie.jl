# In this file lie the function for Gillepsie algorithm


"""
    function give_birth(a::Agent,t,p::Dict)
Used for Gillepsie setting
"""
function give_birth(a::Agent,t,p::Dict)
    new_a = copy(a)
    increment_x!(new_a,t,p)
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
    world[idx_offspring].d = sum(α.(get_x.(skipmissing(world)),Ref(x_offspring))) - α(x_offspring,x_offspring)
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
    function updateWorld_G!(world,t,p)
Updating rule for gillepsie setting.
Returning dt drawn from an exponential distribution with parameter the total rates of events.
 # Args
 t is for now not used but might be used for dynamic landscape
"""
function updateWorld_G!(world,p,update_rates!,t)
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
            world[idx_offspring] = give_birth(mum,t+dt,p)
            update_afterbirth_std!(world,idx_offspring,p)
        end
        return dt
    else
        return -1.
    end
end

"""
    clean_world(world::Array{T}) where {T <: Agent}
Get rid of missing value in `world`
"""
#TODO : put some type specs here
clean_world(world) = collect(skipmissing(world))

"""
    function new_world_G(nagents::Int,p::Dict; spread = 1., offset = 0.)
Returns an array `world0` of size `p["NMax"]`, with `nagents` ancestors agents.
Those agents have traits distributed according to the normal distribution with mean `offset` and standard deviation `spread`.
The dimension of the domain is determined by the size of the array `spread`.
"""
function new_world_G(nagents::Int,p::Dict; spread = 1., offset = 0.)
    typeof(spread) <: Array ? spread = spread[:] : nothing;
    typeof(offset) <: Array ? offset = offset[:] : nothing;
    agent0 = [Agent( spread  .* randn(length(spread)) .+ offset) for i in 1:nagents]
    world0 = vcat(agent0[:],repeat([missing],Int(p["NMax"] - nagents)))
    return world0
end
