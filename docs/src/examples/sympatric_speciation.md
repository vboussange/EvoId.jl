# Modelling Sympatric speciation

This script aims at reproducing the 1999 article of Doebeli [On The Origin of Species By Sympatric Speciation](http://www.nature.com/articles/22521).

In this article, birth and death functions are defined as gaussian, with respective variance ``\sigma_b`` and ``\sigma_d``. It is shown that when ``\sigma_d < \sigma_b``, speciation in the trait space occurs. This is what we reproduce here.

## Running the world
![]()

## Plotting lineages
A cool feature of ABMEv.jl is its ability to track agents ancestors traits (cf [Agent section](../manual/agent.md))

On can plot it, to get an idea of the coalescence of the population.

![]()

Beautiful, isn't it?

!!! tip "Making sense of trait histories"
    Some metrics are available (cf  [Metrics section](../manual/metrics.md)) that summarize the divergence in trait value (or geographical position) through time).
