# In this file lie the function for WF algorithm
"""
$(TYPEDEF)

Wright Fisher algorithm.

In the Wright Fisher process the number of agents is constant through time. It is helpful to visualize it through marbles in jar
![alt text](https://upload.wikimedia.org/wikipedia/commons/0/0b/Random_sampling_genetic_drift.svg)
At each time step, ``N`` agents are picked up from previous generation to reproduce. 
Their number of offspring is proportional to their fitness, calculated as usual with birth and death rates.
It takes thus only one time step to go trough one generation. Thus it is more suitable for numerical simulations than `CFM` or `Gillespie`. 
"""
struct WF <: AbstractAlg end
export WF

function updateWorld_WF!(world,newworld,p,update_rates!,t)
    @debug "updating rates"
    update_rates!(world,p,t);
    # normalise to make sure that we have a probability vector
    fitness = get_fitness.(world)
    # we need to substract the minimum, otherwise fitness can be negative
    prob = normalize(fitness .- minimum(fitness) .+ eps(1.),1)
    offsprings = rand(Multinomial(length(world),prob))
    # this guy can be done in parallel
    @debug "filling new offspring row"
    for (parentID,nboffspring) in collect(enumerate(offsprings))
        for j in 1:nboffspring
            newworld[sum(offsprings[1:(parentID-1)]) + j] = copy(world[parentID])
        end
    end
    # we introduce randomness here
    @debug "incrementing offsprings traits"
    for w in newworld
        increment_x!(w,p)
    end
end
