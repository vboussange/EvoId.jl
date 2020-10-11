# CFM algorithm

Champagnat Ferriere Méléard algorithm described in [Champagnat and Ferriere founding article](https://linkinghub.elsevier.com/retrieve/pii/S0040580905001632). This algorithm similar to Gillepsie algorithms, excepts that it runs much faster!

Indeed, at every time step, only the fitness of the individual picked at random is evaluated. Thus this algorithm is of order ``K`` times more efficient.

In order to use it, you need to feed to the dictionary parameters `p` a constant `Cbar<:Real` that is the upperbound of the maximum of the sum of the birth and death rates (cf article).

!!! warning "Development"
    CFM gives an approximate time step. As of now, we do not manage to obtain qualitatively the same results as the Gillepsie algorithm.
