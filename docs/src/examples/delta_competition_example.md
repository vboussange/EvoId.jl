# A first model of meta population

In this script, we model the evolution of a population where agents are simply defined by their position on some landscape. We implement the simplest possible birth and death function.

## The landscape
Let's start by a linear landscape. We define a discrete segment of length `9`, with reflecting boundary conditions.
