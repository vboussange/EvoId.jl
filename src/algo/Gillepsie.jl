# In this file lie the function for Gillepsie algorithm

"""
$(TYPEDEF)
Gillespie algorithm.

Denote by ``b_i`` and ``d_i`` the birth and death rates of 
agent ``i``. The total rate is given by the sum of all individual rates
``R(t) = \\left[ \\sum_i b_i(t) + d_i(t) \\right]``.
A particular event, birth or death, 
is chosen at random with a probability equal to the rate of this event divided by the total rate ``R``.


# The original article by Gillsepie:
[**A general method for numerically simulating the stochastic 
time evolution of coupled chemical reactions**](https://www.sciencedirect.com/science/article/pii/0021999176900413?via%3Dihub)

"""
struct Gillepsie <: AbstractAlg end
export Gillepsie

function updateBirthEvent!(w::World,::Gillepsie,mum_idx::Int,b,d)
    # updating competition only the two columns corresponding to agent idx
    offspring = give_birth(mum_idx,w)
    x_offspring = get_x(offspring)
    _agents = agents(w)
    for a in _agents
        # Warning : symmetric competition
        tmp = d(get_x(a),x_offspring,w.t)
        a.d += tmp
        offspring.d += tmp
    end
    # Now updating new agent
    offspring.b = b(x_offspring,w.t)
    addAgent!(w,offspring)
end

function updateDeathEvent!(world::World,::Gillepsie,i_event::Int,d)
    x_death = get_x(world[i_event])
    # updating death rate only the two columns corresponding to agent idx
    removeAgent!(world,i_event)
    for a in agents(world)
        a.d -= d(get_x(a),x_death,world.t)
    end
end

"""
    update_rates!(w::World,::Gillepsie,b,d)
This standard updates takes
    - competition kernels of the form α(x,y) and
    - carrying capacity of the form K(x)
"""
function  update_rates!(w::World,::Gillepsie,b,d)
    _agents = agents(w)
    traits = get_x.(_agents)
    # traits = get_xhist.(world)
    n = size(w)
    D = zeros(n)
    # Here you should do a shared array to compute in parallel
    for i in 1:(n-1)
        for j in i+1:n
            C = d(traits[i],traits[j],w.t)
            D[i] += C
            D[j] += C
        end
    end
    # Here we can do  it in parallel as well
    for (i,a) in enumerate(_agents)
        a.d = D[i]
        a.b = b(traits[i],w.t)
    end
end

"""
    function updateWorld!(w::World{A,S,T},g::G,b,d)
Updating rule for gillepsie setting.
Returning `dt` drawn from an exponential distribution with parameter the total rates of events.
"""
function updateWorld!(w::World{A,S,T},g::G,b,d) where {A,S,T,G <: Gillepsie}
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
            updateDeathEvent!(w,g,i_event,d)
        else
            # birth event
            # i_event - n is also the index of the individual to give birth in the world_alive
            updateBirthEvent!(w,g,i_event-n,b,d)
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
