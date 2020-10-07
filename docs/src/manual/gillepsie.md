# Gillepsie algorithm

## Mathematical foundations

- The original article by Gillepsie:
[**A general method for numerically simulating the stochastic time evolution of coupled chemical reactions**](https://www.sciencedirect.com/science/article/pii/0021999176900413?via%3Dihub)


### Rates
Each individual is assigned a birth $` b_i `$ and death $` d_i `$ rate. The total rate is given by the sum of all individual rates
```math
R(t) = \left[ \sum_i b_i(t) + d_i(t) \right]
```
A particular event, birth or death, is chosen at random with a probability equal to the rate of this event divided by the total rate $`R`$
> This has to be checked, we are still not hundred percent sure

### Time steps
An event is exponentiallly distributed in time, with parameter $`\lambda = U(t)`$. This makes events memoryless, meaning that the probability of having a birth or death event is always the same, no matter when ($`P(X > s_t | X > t) = P(X > s) `$.
> Let $`B(t) = \sum_i b_i(t)`$ and  $`D(t) = \sum_i d_i(t)`$. Let $`T_b, T_d`$ the time for a birth or death event to occur. Then we have $`P(T_b < T_d) = \frac{B(t)}{B(t) + D(t)}`$  (competing exponentials).
#### Inversion method
Let $`U`$ be an $`\mathcal{U}_{(0,1)}`$-distributed random variable and $`F \colon \R \to [0,1]`$ be a distribution function. Then we have
```math
P(I_F(U) \leq x ) = P(U \leq F(x)) = F(x)
```


Thanks to the ***inversion method*** we get the incremental time step $`dt`$, exponentially distributed with parameter $`\lambda = R(t)`$, as

```math
    dt(\omega) = -\frac{\log(U(\omega))}{R(t)} \iff X(\omega) = \exp(-U(t)dt(\omega))
```
## Initialize
```@docs
new_world_G
```
## Run
```@docs
runWorld_store_G
```

## Scenarios
As of now, no mode is implemented. For further examples, check the folder `examples` in source code.

## Specific parameters

- ```dt_saving = 10.```
will allow to save the world every 10. time steps. If not specified, the algorithm will return first and last time step world.
- ```NMax``` Maximum number of individuals that can be attained. If attained, then the programm stops.

## Developping
### Efficiency
The simulation are still very long.
:flushed: How to improve it?
- We think it would be more efficient if we found an other way of incrementing mutations


### Parallelism
> For now there is no parallelism implemented for one run

> But we think we should rather set up a pmap or the macro `@Threads` to explore parameter space


```@autodocs
Modules = [ABMEv]
Pages   = ["ABMEv_Gillepsie.jl"]
```
