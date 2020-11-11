# Taken from gillepsie algorithm documentation

!!! tip "Inversion method"
    Let ``B(t) = \sum_i b_i(t)`` and  ``D(t) = \sum_i d_i(t)``. Let ``T_b, T_d`` the time for a birth or death event to occur. Then we have ``P(T_b < T_d) = \frac{B(t)}{B(t) + D(t)}``  (competing exponentials).

    Let ``U`` be an ``\mathcal{U}_{(0,1)}``-distributed random variable and ``F \colon \R \to [0,1]`` be a distribution function. Then we have
    ```math
    P(I_F(U) \leq x ) = P(U \leq F(x)) = F(x)
    ```


    Thanks to the ***inversion method*** we get the incremental time step ``dt``, exponentially distributed with parameter ``\lambda = R(t)``, as

    ```math
        dt(\omega) = -\frac{\log(U(\omega))}{R(t)} \iff X(\omega) = \exp(-U(t)dt(\omega))
    ```

    ### Time steps
    An event is exponentiallly distributed in time, with parameter ``\lambda = U(t)``. This makes events memoryless, meaning that the probability of having a birth or death event is always the same, no matter when (``P(X > s_t | X > t) = P(X > s) ``.
