# Gillepsie simulation

## Diversification
### Gaussian birth coefficient, Constant carrying capacity
```julia
a = 0;
sigma_K = .9;
sigma_a = .7;
K0 = 1000;
K(X) = gaussian(X[1],0.,sigma_K)
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
# α(X,Y) = 0.
p_default = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 2000,
        "NMax" => Int(2000),
        "dt_saving" => 10.0)
na_init = K0
world0 = new_world_G(na_init,p_default,spread = .01)
```
![EVOID_div_t4000](uploads/3a49ff4fe4db161bf360eea97694ff26/EVOID_div_t4000.png)
### Constant birth coefficient, Gaussian carrying capacity
```julia
a = 0;
sigma_K = .9;
sigma_a = .7;
K0 = 1000;
K(X) = 1
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0 / gaussian(X[1],0.,sigma_K)
# α(X,Y) = 0.
p_default = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 4000,
        "NMax" => Int(2000),
        "dt_saving" => 20.0)
na_init = K0
```
![EVOID_bis_div.ong](uploads/8e1f821923afd74902b3ec6567a1736d/EVOID_bis_div.ong.png)
> what you could do would be to plot the adaptive dynamics of the monomorphic populations

### Quadratic birth rate
```julia
a = 0.125;
sigma_a = 1.251;
K0 = 1000;
K(X) =  1 - a * X[1]^2
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
# α(X,Y) = 0.
p_default = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 4000,
        "NMax" => Int(2000),
        "dt_saving" => 20.0)
na_init = K0
```
![Unknown](uploads/e6b6bc135fbd2f48cc5ad8bd18854420/Unknown.png)

#### Equivalence

![PDE_quad_termsol](uploads/c1f5c670a9d82df1349ed473b9954135/PDE_quad_termsol.png)
![EVOID_quad_time_average_distrib_deep_time](uploads/2613b1f5e919fdbc7ee3386c680a1908/EVOID_quad_time_average_distrib_deep_time.png)

## No diversification
```Julia
a = 0;
sigma_K = .9;
sigma_a = 1.0;
K0 = 1000;
K(X) = gaussian(X[1],0.,sigma_K)
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
# α(X,Y) = 0.
p_default = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 2000,
        "NMax" => Int(2000),
        "dt_saving" => 10.0)
na_init = K0
```
![Gillepsie_results_nodiversification](uploads/3a459012508c85cf853246a37537f160/Gillepsie_results_nodiversification.png)
