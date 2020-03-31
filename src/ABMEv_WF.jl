# In this file lie the function for WF algorithm
    """
        function updateWorld_WF!(world,newworld,C,p,update_rates!,t,reflected)
            If reflected=true, we reflect only first trait corresponding to geographic position
    """
    function updateWorld_WF!(world,newworld,p,update_rates!,t,reflected)
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
            increment_x!(w,p,reflected=reflected)
        end
    end
