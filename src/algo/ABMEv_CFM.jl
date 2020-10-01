# CFM = Champagnat Ferriere Meleard
# This contains all methods for algorithm described in article
# DOI : 10.1016/j.tpb.2005.10.004

# In this file lie the function for Gillepsie algorithm
"""
$(TYPEDEF)
"""
struct CFM <: AbstractAlg end
export CFM

function updateBirthEvent!(w::World,::CFM,mum_idx::Int)
    # updating competition only the two columns corresponding to agent idx
    offspring = give_birth(mum_idx,w)
    addAgent!(w,offspring)
end

function updateDeathEvent!(world::World,::CFM,i_event::Int)
    removeAgent!(world,i_event)
end

"""
$(SIGNATURES)
Updating rule for CFM setting.
algorithm described in article
DOI : 10.1016/j.tpb.2005.10.004
"""
function updateWorld!(w::World{A,S,T},c::CFM) where {A,S,T}
    # total update world
    alive = agents(w)
    events_weights = ProbabilityWeights(vcat(get_d.(alive),get_b.(alive)))
    # Total rate of events
    ∑ = sum(events_weights)
    dt = - log(rand(T))/∑
    update_clock!(w,dt)
    if dt > 0.
        i_event = sample(events_weights)
        # This is ok since size is a multiple of 2
        n = size(w)
        if i_event <= n
            # DEATH EVENT
            # In this case i_event is also the index of the individual to die in the world_alive
            updateDeathEvent!(w,g,i_event)
        else
            # birth event
            # i_event - n is also the index of the individual to give birth in the world_alive
            updateBirthEvent!(w,g,i_event-n)
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
