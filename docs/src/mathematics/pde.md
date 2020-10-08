# PDE Equivalence
# Equivalence with PDE models
> HDR Champagnat Theorem 2.5 Section 2.2

> Recall that the IBM is a Markov process in $`\mathcal{M}^1`$ which state at time $`t`$ is the measure
```math
    \nu_t = \sum_{i=1}^{N_t} \delta_{x_i}
```

> Hence an individual with trait $`x`$ in population $`\nu_t`$ gives birth to a new individual at rate $`b(x)`$ and dies at rate
```math
    d(x) + \int c(x,y)(\nu_t(dy) - \delta_x(dy)) = d(x) - c(x,x) + \sum_{i=1}^{N_t} c(c,x_i)
```
The Dirac mass in the integral stands for the fact that an individual is not in competition with itself. Hence when $`N_t = 1`$ the competition term cancels. When there is birth,  with probability $`p(x)`$ the offspring has trait $`y = x + H`$ where $`H`$ is a random variablee with law $`m(x,h)\,dh`$.

## Assumptions
Assume that
#### Assumption 1 : bounds and parameters regularity.
- Functions $`b,c,d,p`$ are continuous and bounded, with positive are nul values
- $`\exists \bar{m}, \forall x \in D, h \in \R^d m(x,h) \leq \bar{m}(h)`$

#### Assumption 2
$`\nu_t^K = \frac{1}{K}\nu_t`$ where $`\nu_t`$ is constructed such that $`c\equiv \frac{1}{K}c`$


 

Then we get the **large population limit without mutation scaling**

```math
    \partial_t u(t,x) = u(t,x) \left((1-p(x))b(x) -d(x) \int_D c(x,y) u(t,y) dy\right) + \int_D b(x) \, p(y) \, u(t,y) \, m(y,x-y)dy
```


## Mutations
Assuming a small mutational variance $`\max \sigma^2_i << 1`$ and a mutation rate $`U`$, the mutational effects can be approximated by an elliptic operator $`\sum_{i=1}^{n} (\mu_i^2/x)\partial_{ii}`$ with $`\mu_i = \sigma_i\sqrt{U}`$
> :warning: check that with Burger

In other words (from Champagnat, Ferriere and Meleard 2006), we have

```math
\partial_t u(t,x) =  u(t,x)\big(b(t,x) - d(x,u(t,x)) \big) + \frac{1}{2} \Delta (\sigma^2 r \mu u)(t,x)
```
where 
```math
d(x,u(t,x)) = (u \ast c)(x)
```
and ``c`` is the competition function.
For us, ``x \in ``