# Callbacks

Callbacks can be used to extract properties of the world at each `dt_saving` time steps of your simulation.


## Constructing the callbacks
A callback has to be of the form

```julia
cb = (names = String[], agg = Function[])
```
It is a tuple, with first value corresponding to the names of the aggregate properties of the world. The second correspond to the aggregation functions.

We provide here an example on how to extract the ``\gamma`` diversity of a simulation biological population. ``\gamma`` diversity can be calculated as the variance of the trait distribution of the population.
Here is how we write the function
```julia
cb = (names = ["gamma_div"], agg = Function[w -> var((get_x(w,1)))])
```

Here is how we use it
```julia
myspace = (RealSpace{1,Float64}(),)
sigma_K = .9;
sigma_a = .7;
K0 = 1000;
b(X) = gaussian(X[1],0.,sigma_K)
d(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
D = (1e-2,)
mu = [.1]
NMax = 10000
tend = 1.5
p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax

myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
w1 = copy(w0)
@time sim = run!(w1,Gillepsie(),tend,cb=cb,dt_saving = .1)
```

## Accessing the callbacks
You can easily access the properties, using `sim` as you would for a usual `Dictionary`.

```julia
using Plots
plot(get_tspan(sim),sim["gamma_div"])
```
