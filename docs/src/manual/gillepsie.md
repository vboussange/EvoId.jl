# Gillepsie algorithm

## Mathematical foundations

- The original article by Gillepsie:
[**A general method for numerically simulating the stochastic time evolution of coupled chemical reactions**](https://www.sciencedirect.com/science/article/pii/0021999176900413?via%3Dihub)


### Update Rates
 `` b_i `` and `` d_i `` represent respcetively birth and death rates of agents ``i``. The total rate is given by the sum of all individual rates
```math
R(t) = \left[ \sum_i b_i(t) + d_i(t) \right]
```
A particular event, birth or death, is chosen at random with a probability equal to the rate of this event divided by the total rate ``R``.


```@autodocs
Modules = [EvoId]
Pages   = ["EvoId_Gillepsie.jl"]
```
