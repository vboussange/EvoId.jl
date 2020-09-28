# In this file lie the function for Gillepsie algorithm

struct Gillepsie <: AbstractAlg end

"""
    function give_birth(a::Agent,t,p::Dict)
Used for Gillepsie setting
"""
function give_birth(mum_idx::Int,w::World{AbstractAgent{A,R},S,T}) where {A <: Ancestors{true},R,S,T}
    new_a = copy(w[mum_idx])
    increment_x!(new_a,space(w),parameters(w),time(w))
    return new_a
end

function updateBirthEvent!(world::World,::Gillepsie,mum_idx::Int)
    # updating competition only the two columns corresponding to agent idx
    @unpack d,b = parameters(world)
    offspring = give_birth(w,mum_idx)
    x_offspring = get_x(offspring)
    _agents = agents(world)
    for a in _agents
        a.d += d(get_x(a),x_offspring)
    end
    # Now updating new agent
    offspring.d = sum(d.(get_x.(_agents),Ref(x_offspring))) - d(x_offspring,x_offspring)
    offspring.b = b(x_offspring)
    addAgent!(w,offspring)
end

function updateDeathEvent!(world::World,::Gillepsie,i_event::Int)
    @unpack d,b = parameters(world)
    x_death = get_x(world[i_event])
    # updating death rate only the two columns corresponding to agent idx
    removeAgent!(world,i_event)
    for a in agents(world)
        a.d -= d(get_x(a),x_death)
    end
end

"""
    function updateWorld_G!(world,t,p)
Updating rule for gillepsie setting.
Returning dt drawn from an exponential distribution with parameter the total rates of events.
 # Args
 t is for now not used but might be used for dynamic landscape
"""
function updateWorld!(w::World,g::G) where {G <: Gillepsie}
    # total update world
    agents = agents(world)
    events_weights = ProbabilityWeights(vcat(get_d.(agents),get_b.(agents)))
    # Total rate of events
    U = sum(events_weights)
    dt = - log(rand(T))/U
    update_clock!(world,dt)
    if dt > 0.
        i_event = sample(events_weights)
        # This is ok since size is a multiple of 2
        n = size(world)
        if i_event <= n
            # DEATH EVENT
            # In this case i_event is also the index of the individual to die in the world_alive
            updateDeathEvent!(world,G,i_event)
        else
            # birth event
            # i_event - n is also the index of the individual to give birth in the world_alive
            mum = world[i_event-n]
            update_afterbirth_std!(world,idx_offspring)
        end
        return dt
    else
        return -1
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
