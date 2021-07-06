# Plotting

EvoId comes with Plot recipes.
```julia
function plot(sim::Simulation;trait = 1)
```

## Arguments
- if `length(trait) == 1` then we scatter plot `trait` along time
- if `2 <= length(trait) <= 3` then we project world of the last time step in the two  (three) dimensional trait space define by `trait`

!!! warning "To be implemented"
        We might want to get a 3 dimensional scatter plot
        with time, trait1 and trait2 as axis
