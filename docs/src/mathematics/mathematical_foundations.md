# Mathematical foundations

- [Dieckmann and Law](https://link.springer.com/article/10.1007%2FBF02409751) (probably more accessible than following), and then
- [Nicolas, Ferriere and Meleard](https://www.sciencedirect.com/science/article/pii/S0040580905001632?via%3Dihub)

showed that the canonical equations of adaptive dynamics, describing the evolution in the phenotypic space, can be derived by considering the stochastic individual-based model corresponding to
```math
\partial_t u(t,x) = u(t,x)(1-\frac{\int \alpha(x,y) u(y,t) dy}{K(x)})
\tag{1}
```
in the limit of **rare mutations, small mutational effects, and infinite population sizes**.

Under these assumptions, Dieckmann and Law showed that adaptive dynamics is the first-order approximation of the mean path averaged over infinitely many realizations of the stochastic simulations obtained from the individual-based model.

[[_TOC_]]

## The basics of Adaptive Dynamics
> This section is inspired from *Adaptive Diversification, Doebeli 2011*

Imagine a monomorophic population with trait ``x``, which follows the dynamics
```math
\partial_t u(t,x) = u(t,x)\big[b(x) - c(x)u(t,x)\big]
```
> Logistic models are thought to be mathematically representative of a large class of models (During2008) :smirk:. Doebeli "Chaos" paper might also be relevant for this purpose


Equilibrium is given by ``u^*(x) \equiv K(x) = b(x)/c(x)``. Hence one can reformulate equation above by
```math
    \partial_t u(t,x) = u(t,x) \left(b(x) - \frac{b(x)u(t,x)}{K(x)}\right )
```

Adaptive dynamics aims at determining the evolutionary trajectory of the trait ``x`` by considering the fate of rare mutants with traits ``y`` in resident monomorphic population with trait ``x``. The resident population is assumed to be at equlibrium ``K(x)``. Because mutant is rare, the mutant's population dynamics is only affected by the density of the resident, which in turn is unaffected by the mutant's invasion attempt, and hence remains at ``K(x)``. Hence the effective density that mutant experiences during invasion attempt is detrmined by ``\alpha(x,y)K(x)``. The function  ``\alpha`` is the **competition kernel** and describes the strength of competition that exerts the phenotype ``x`` on phenotype ``y``.
> For our purpose, we consider it as symmetric, which corresponds to the canonical biological examples of birds with different beak size. We think this is also valid for plants, and in general along a genus.

> It is assumed that ``\alpha(x,x) = 1`` (scaled to unity).

#### Stabilising selection
Stabilising selection ensures that trait ``x`` is contained, thus avoiding regions of extreme trait values that would be biologically unrealistic. It can appear by limiting birth rate ``b(x)`` or increasing death rate ``c(x)`` for extreme x

Now we can determine the population dynamics of a mutant with population size ``N_{mut}`` and trait ``y`` experiencing the effective density ``N_{eff} = \alpha(x,y)K(x)``
```math
\begin{aligned}
    \frac{\partial N_{mut}}{\partial x}(t,x) &= N_{mut}(b(y) - c(y)N_{eff})\\
    &= N_{mut}(b(y) - c(y)\, \alpha(x,y)K(x))\\
    &= N_{mut}(b(y) - \frac{b(y)\, \alpha(x,y)K(x)}{K(y)})
\end{aligned}
```

## Invasion fitness
The invasion fitness `` f(x,y) `` corresponds to the fitness of a mutant with traits ``y`` in a resident monomorphic population with trait ``x``
```math
\begin{aligned}
    f(x,y) &= \frac{1}{N_{mut}}\frac{\partial N_{mut}}{\partial t}  = \left( b(y) - c(y)N_{eff} \right)\\
    \iff \\
    f(x,y) &= b(y)\left( 1 - \frac{\alpha(x,y)K(x)}{K(y)} \right)
\end{aligned}
```
> :gun: Plot the fitness function in the case of the EVOID
## Canonical equation for the Adaptive Dynamics
```math
\frac{d\, x}{d \, t} = m(x) D(x)  \tag{2}
```
> You should write it in the multidimensional case
- ``D(x) = \frac{\partial f(x,y)}{\partial y}\Big|_{y=x}`` is called the **selection gradient**
- ``m(x)`` is the **mutational kernel**. In Champagnat2011, it is denoted as ``m(x) = \mu(x) \, \frac{\sigma_0^2(x)}{2} \, N(x)`` where ``\mu`` denottes the probability that a birth from an individual with trait ``x`` gives rise to a mutation. ``\sigma_0^2(x)`` denotes the variance of the distribution of a mutant trait ``y`` born from an individual with trait ``x``. Even if the mutational effect is independent of ``x``, the rate at which new mutations occur depends on the current population size ``N(x)``. If population is monomorphic, ``m`` only scales time. However, when population becomes polymorphic, mutation affects evolutionary dynamics.

> :question:Essentially, adaptive dynamics is a first order approximation of the nonlinear dynamics of any evolutionary model
### Symmetric competition kernels
Here we assume that ``b(x) = b \in \R`` and ``c(x) = b/K(x)``. Thus the canonical equation for adaptive dynamics yields
```math
    \frac{d\, x}{d \, t} = b\left( \alpha(x,y)\frac{ K' (x)}{K(x)} - \partial_2 \alpha(x,y) \right).
```
Assuming that ``\alpha(x,y) := \alpha(|x-y|)`` and maximum competition at 0, we have that ``\partial_2 \alpha(x,x) = 0`` and thus
```math
    \frac{d\, x}{d \, t} = b \, \alpha(x,y)\frac{ K' (x)}{K(x)}
```
#### :heart: Evolutionary singularities: ``D(x^*) = \partial_2 f(x^*,x^*) = 0``
The fixed points of equation (2) are the points where the fitness gradient nullifies, and are called *evolutionary singularities*.

> **In this special case, the monomorphic population evolutionary dynamics is only driven by carrying capacity.** Assuming a unimodal carrying capacity, the adaptive dynamics will always converge to the maximum of ``K(x)``

> :question: This can be made explicit using the concept of singular points and convergence stability

>:rocket: unimodality of K(x) implies the existence of a unique singular point at the trait value maximising the carrying capacity
****


#### :heart: Convergence stability: ``D'(x^*) < 0``
Convergence stability of the singular point is determined by the quantity ``D'(x^*)``

> Because ``x^*`` is a maximum, we have ``D'(x^*) = b\frac{K''(x)}{K(x)} < 0`` hence the singular point is always convergence stable.

By definition, selection vanishes at the singular point. Therefore, after convergence, **second order effects of selection come at play**, which is why the competition kernel and hence frequency dependence starts to be important.

 In a ``d`` dimensional setting, convergence stability is obtained if the Jacobian matrix of the set of ODE (2) has all its eigenvalues negative at the singular point ``x^*``.
 > :warning: in a multidimensional space, convergence stability nor being in the neihborhood of a singular point are required for evolutionary branching to occur.
#### :heart: Evolutionary branching point: ``\partial_{22}f(x^*,x^*) < 0``


> In 1D, convergent stable singular point is either a maximum of minimum of the invasion fitness function. Note that for any ``y \neq x^*``, ``\alpha(x^*,y)<\alpha(x^*,x^*)`` and ``K(y)<K(x^*)``.

:rocket: Whether ``x^*`` can be invaded by a mutant depends on the relative magnitude of these two effects, which in turns is determined by the **curvature of the competition kernel and carrying capacity**.

> In our particular case we have
```math
\partial_{22}f(x^*,x^*) = b \left(\frac{K''(x^*) }{K(x)} - \alpha''(x^*,x^*) \right)
```
[]()

> Hence we get that ``x^*`` is a fitness minimum, i.e. **an evolutionary branching point**``\iff``
```math
\alpha''(x^*,x^*)<\frac{K''(x^*) }{K(x)}
```
##### Reminder:
> ![Alt Text](img/maxmin2.gif)

:exclamation: in dimension D, the condition for minima is obtained by the determinant of the Hessian ([The second partial derivative test](https://en.wikipedia.org/wiki/Second_partial_derivative_test))

> Note that if  ``K(x) = \exp(-\frac{(x-x_0)^2}{2\sigma_K^2})`` and ``\alpha(x,y) = \exp(-\frac{(x-y)^2}{2\sigma_\alpha^2})``
then above condition boils down to
```math
    \sigma_\alpha < \sigma_K
```

>The Appendix of Doebeli2011 demonstrates from general dynamic systems theory that phenotypes ``x_1`` and ``x_2`` can coexist if they are close enough to the singular point.
#### Taylor expansion of the fitness function
##### 1D
Let's expand the fitness function for a scalar trait
```math
    f(x,y) = f(x,x) + \partial_2 f(x,x)(y-x) + \partial_{22} f(x,x)\frac{(y-x)^2}{2} + \dots
```
First term is 0 for all ``x``. The second term vanishes at the singular point, which is where the second order term come at play. If ``\partial_{22} f(x^*,x^*) <0`` no nearby mutants can invade the resident population that is monomorphic and thus conditions for evolutionary stability are obtained.

In contrast, ``\partial_{22}f(x^*,x^*) > 0 `` is the condition for evolutionary instability, or potential evo- lutionary branching points, as the mutant now can potentially invade the resident.

##### 2D
```math
    f(x,y) = f(x,x) + \nabla_2 f(x,x)(y-x) + \frac{1}{2}(y-x)^{T} \nabla_{2}^2 f(x,x) (y-x) + \dots
```
We refer to ``\nabla_{2}^2 f(x,x) \equiv H(x)`` as the hessian matrix.

In this scenario, the first-order term can become zero or be
arbitrarily close to zero when the vector ``y-z`` is orthogonal or nearly orthogonal to the gradient ``\nabla_2 f(x,x)``. In other words, the second-order terms for trait values y that lie orthogonal to the direction of the gradient of f(x, y) become significant regardless of whether or not the trajectory is in the vicinity of a singular point. In particular, if H(x)is negative definite, no nearby mutants can invade the resident population that is monomorphic and we have the conditions for evolutionary stability.

:exclamation: When the Hessian matrix has an eigenvector with a positive eigenvalue and this eigenvector has a nonzero projection in the subspace orthogonal to the direction of the selection gradient, we get the condition for evolutionary instability.
The mutants from these directions can potentially invade the resident. Geometrically, the Hessian matrix at a particular point determines the local curvature of the level set of the function ``f(x, y)`` as a function of y passing through that point. Depending on this curvature (the non- negativity of the Hessian matrix), there may exist directions along which the invasion fitness function has a minimum with respect to the mutant trait values and hence becomes evolutionary unstable.

#### Determining the dimorphic steady state: ``D({\bf x}) = D(x_1,x_2) = \begin{pmatrix}D_1(x_1,x_2) \\ D_2(x_1,x_2)\end{pmatrix} = \begin{pmatrix}\partial_3f(x_1,x_2,x_1)\\\partial_3f(x_1,x_2,x_2) \end{pmatrix} = {\bf 0}``
Consider mutants in the dimorphic population ``(x_1,x_2)``. Ecological dynamics yields
```math
    \begin{aligned}
        \partial_t N_1 = bN_1(1-\frac{N_1 + \alpha(x_1,x_2)N_2}{K(x_1)}) \\
        \partial_t N_2 = bN_2(1-\frac{N_2 + \alpha(x_1,x_2)N_1}{K(x_2)})
    \end{aligned}
```
Equilibrium densities are given by
```math
    \begin{aligned}
    N_1^* = N_2^* = \frac{\exp\left(x_1^2\left( \frac{2}{\sigma_\alpha^2} - \frac{1}{2\sigma_K^2} \right)\right)}{1+\exp\left(\frac{2x_1^2}{\sigma_\alpha^2}\right)}
    \end{aligned}
```
Hence invasion fitness in the dimorphic population ``(x_1,x_2)`` is
```math
    f(x_1,x_2,y) =b\left(1 - \frac{\left( \alpha(x_1,y)N_1^* + \alpha(x_2,y)N_2^*\right)}{K(y)}\right)
```
Thus the selection gradient in the resident ``x_1`` yields
```math
    \begin{aligned}
    D_1(x_1,x_2) &= \partial_3f(x_1,x_2,x_1)\\
    &= \dots
    \end{aligned}
```
The 2-dimensional adaptive dynamics systems yields
```math
    \frac{d {\bf x}}{dt} = m({\bf x})D({\bf x})
```
where the mapping ``m \colon \R^2 \to\R^2`` decsribes the mutational process in the two resident strains.
> look Appendix Doebeli2011 that deals with mutations

With symmetric Gaussian competition we get
```math
    x_1^* = x_2^* = \sqrt{\ln\left(\frac{2\sigma_K^2}{\sigma_\alpha^2} - 1\right)\frac{\sigma_\alpha^2}{2}}
```
> :gun: This is a good way to check that your IBM is working!!!

#### Convergence stability of the dimorphic state: ``J_D(x_1^*,x_1^*)``
If the eigenvalues of the Jacobian matrix ``J_D(x_1^*,x_1^*)``are negative, then after after evolutionary branching the two phenotypic branches converge to ``x_1^*`` and ``x_2^*``.

#### Evolutionary branching of the dimorphic state
:flushed: because of the symmetry assumption the two singular strategies either both represent local maxima or minima for the invasion fitnessfunction. It can be shown that whenever ``\sigma_\alpha<\sigma_K`` both singular strategies are fitness minima, hence **evolytionarily unstable**.

#### Behaviour in deep time
Evolutionary branching continues ad infinitum, that is when time goes to inifinity there is infinitely many branches with decreasing population starting from carrying capacity maximum. This is a peculariaty of the Gaussian kernel. Hence the popupulation will simply be polymorphic, normally distributed. This is what we observe from our PDE model.

> :exclamation: Doebeli2011 claims that Gaussian ecological functions in stochastic IBM for finite populations typically leads to only a small number of successive events

> :exclamation: This is why it makes sense to extend this work to a continuous setting.

### Quadratic carrying capacity
Competition kernel is still assumed gaussian. Now
```math
    K(x) = 1-a\,x^2
```
Conditions for the singular point to be an evolutionary branching point yields
```math
    \sigma_\alpha < \frac{1}{\sqrt{2a}}.
```
Not possible to get explicit expression for the singular coalition of the two dimensional dynamics, but one can give examples and solve it numerically. The coexisting strategies are evolutionary stable as long as ``\sigma_\alpha`` is not too small.
For example, take
```
a= 0.125
```
Evolutionary stability is obtained for ``1.25 < \sigma_\alpha < 2 ``
> :flushed: also important to note that even if the coexisting points are not convergent stable the four dimensional adaptive dynamics can have a singular coalition in which all strains are evolutionary stable.

:exclamation: **This proves that even if frequency dependent is strong enough to induce diversification, it noes not necessarily lead to an infinite series of subsequent branching events.**

>:gun: We made the simulation for the corresponding PDE, and this is what we get
> ![Alt Text](img/adaptive_rad_1.png)
> Note that is proves to be a stiff problem when ``\mu \neq 0``

#### Derive the particular case for Eq. (1)

### Particular case of the logistic map
At equilibrium, the monomorphic resident population is distributed along the resource such that  ``u(t,x) = K(x)`` since we consider that ``\alpha(x,x)`` = 1.
> - This is different for mutator selector equation, but you should check it
> - Also check it with varying birth rate, as Champagnat example
Invasion fitness for the logistic map, of a rare mutant ``y`` is its per capital rate of growth in a resident population with phenotype ``x`` and is given by
```math
 f(x,y) = \left(1-\frac{\alpha(x,y)K(x)}{K(y)} \right)
```
Indeed
> This has to be checked

## Competition and resource Kernel
- [**Evolutionary dynamics from deterministic microscopic ecological processes**](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.101.032411) shows that adding a power to the kernels does not modify the adaptive dynamics evolutionary trajectories. It just accelerates its rates.
- [**The shape of the competition and carrying capacity kernels affects the likelihood of disruptive selection**](https://www.sciencedirect.com/science/article/pii/S0022519309000733?via%3Dihub) presents alternative functional forms for competition and resource kernels, and investigate its impact on diversification. Box-like kernel can facilitate evolutionary branching.

## Conditions for diversification
Under the dynamics above, we have that adaptive diversification happen under the condition
```math
\frac{\partial^2 f }{\partial x^2}(x,y)
```
> This is also tackled in Champagnat's memoire is Section 1.4.3, with example in Section 1.2.2.

## Mutation
Here we discuss how we can introduce mutation in Eq. (1) so that we obtain

```math
\partial_t u(t,x) = u(t,x)(1-\frac{\int \alpha(x,y) u(y,t) dy}{K(x)}) + \Delta_x u(t,x)
\tag{mutations}
```



## Correspondance with Individual  Based Model
[The Canonical Equation of Adaptive Dynamics: A Mathematical View](http://www.akademiai.com/openurl.asp?genre=article&id=doi:10.1556/Select.2.2001.1-2.6) has derived for the first time canonical expression from first principles, based on taking limits of a jump process
- inifinite population
- mutation are very rare, that is **ecological and evolutionary timescales** are separated.
### References
- A very good introduction in Adaptive Dynamics on Wikipedia (in French):
[**Dynamique Adaptive**](https://fr.wikipedia.org/wiki/Dynamique_adaptative) 
