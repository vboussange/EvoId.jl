# The geotrait
> This note has been taken from `2Y/articles/geotrait/mathematical_notes.md`

!!! note "'note' admonition"
    In Landscape Ecology, it is of particular interest to study the geographical position of the lineage through time. This is because richness patterns are thought to arise through allopatric speciation, where populations get separated in space and time. This is the topic of the following section.

## Julia accessors
An agents of type `Agent{Ancestors{true}}` stores the values of its ancestors traits. From it, one can thus study the history of the lineage.

## Projection of the geographic history

Consider a set of individuals which lineage geographical history has been stored in a vector ``x^{(i)}(t)``. If individual ``i`` was born at time ``t^*`` then ``x^{(i)}(t>t^*)`` represents its geographical position, while ``x^{(i)}(t<t^*)`` represents its ancestors geographical position.

### Isolation in time
We want to account for the isolation in time and possibly in space of the lineage of a given individual. That is, for how long and how distant stay lineages appart?

#### Setting
Consider a discrete setting, where ``N(t)`` individuals (or populations) evolve over a set of ``M`` demes disposed in a linear fashion, such that ``x^i(t) \in \{1,2,...,M\}``. The population is characterised by the counting measure
```math
\nu = \sum_i^{N(t)} \delta_{x^i(t)}
```


We define the geographical history hamming distance ``\mathfrak{g}`` as
```math
\mathfrak{h}\big(x^{(i)}(t),x^{(j)}(t)\big) = \int_0^t \text{ceil}(\frac{|x^{(i)} - x^{(j)}|(s)}{M-1})ds
```

We extend this definition with the measure ``\mathfrak{h}^*`` which takes into account geographic distance
```math
\mathfrak{h}^*\big(x^{(i)}(t),x^{(j)}(t)\big) = \int_0^t \Big[x^{(i)} - x^{(j)} \Big]^2(s)\, ds
```

Finally, we introduce the geotrait distance ``\mathfrak{g}`` as
```math
\mathfrak{g}\big(x^{(i)}(t),x^{(j)}(t)\big) = \Big[\int_0^t (  x^{(i)} - x^{(j)} )(s) \, ds\Big]^2
```

Note that by the triangle inequality we have that

```math
\mathfrak{g}\big(x^{(i)}(t),x^{(j)}(t)\big) \leq \mathfrak{h}^*\big(x^{(i)}(t),x^{(j)}(t)\big)
```
Equality should arises if ``x^{(i)}, x^{(j)}`` are positively linearly dependent (which should not be the case).

However, what we are eventually interested is a population measure. This measure should be related to the average pairwise distance across the population. Hence we define ``\mathcal{D_d}(\nu,t)`` such that

```math
\mathcal{D_d}(\nu,t) = \frac{1}{2N^2}\sum_i^N  \sum_j^N d(x^{(i)},x^{(j)},t)
```

Let ``g^{(i)}(t) = \int_0^t x^{i}(t) dt`` and ``G(t) = \{g^{(i)}(t), i \in \{1,2,\dots,N(t)\}\}``. By observing the following

```math
\frac{1}{2N^2}\sum_{i,j}^N (y_i-y_j)^2 \\=
 \frac{1}{2N^2}\sum_{i,j}^N ((y_i-\bar{y}) -(y_j-\bar{y}))^2 \\
 = \frac{1}{2N^2} 2N \sum_i^N(y_i-\bar{y})^2 = \text{Var}(Y)
```

Thus we have ``\mathfrak{g}\big(x^{(i)}(t),x^{(j)}(t)\big) = [g^{(i)}(t) - g^{(j)}(t)]^2``  and hence
```math
    \mathcal{D_\mathfrak{g}}(\nu,t) = \text{Var}(G).
```

Also  remark that

```math
\mathcal{D_{h^*}}(\nu,t) = \frac{1}{2N^2}\sum_i^N  \sum_j^N h^*(x^{(i)},x^{(j)},t) \\
= \frac{1}{2N^2}\sum_i^N  \sum_j^N  \int_0^t \Big[x^{(i)} - x^{(j)} \Big]^2(s)\, ds \\
= \int_0^t  \frac{1}{2N^2}\sum_i^N  \sum_j^N  \Big[x^{(i)} - x^{(j)} \Big]^2(s)\, ds \\
= \int_0^t \text{Var}(X)(s)ds

```

One could also imagine a value ``h^{(i)}(t) = \frac{1}{2N}\sum_j^{N(t)}\mathfrak{h}^*\big(x^{(i)}(t),x^{(j)}(t)\big)`` and in this case we would have
```math
    \mathcal{D_\mathfrak{h}}(\nu,t) = \frac{1}{N}\sum_i^{N(t)} h^{(i)}(t)
```

### Mobility
How much do lineages move?

Here is a time average of the speed
```math
    \frac{1}{N} \sum_i < \partial_t x^{(i)}(t) >_t \\
    =  \frac{1}{N} \sum_i \int_0^t \partial_s x^{(i)}(s) ds \\
    = \frac{1}{N} \sum_i [x^{(i)}(t) - x^{(i)}(0))]
```

But one could also have a moving average, that is, averaging
```math
    \frac{1}{N} \sum_i \sum_j \int_{j\tau}^{j(\tau+1)} \frac{1}{\tau} \partial_s x^{(i)}(s) ds
```
