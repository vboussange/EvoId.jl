"""
    run!(w, alg, tend, b, d; dt_saving=nothing, cb=(names = String[],agg =nothing))

Run `w` with algorithm `alg`, until `tend` is reached. User needs to provide `b` the birth function,
which takes as arguments `X,t`, and provide `d` the death function, with arguments `X,Y,t`.
The run stops if the number of agents reaches`p["NMax"]`, where `p` is the parameter dictionary in the world `w`.


Returns a `Simulation` object. It is a container for snapshots of the world at every `dt_saving` 
time steps. It renders post processing easier, through dedicated methods to obtain time series of quantities.

# Keyword arguments
- if `dt_saving` specified, world is saved every time steps.
If `dt_saving` not specified, `sim` contains an array of two elements,
first corresponding to initial conditions and last corresponding to world in the last time step.
- if `t_saving_cb` specified, callbacks are computed at each steps time specified in the array.
This functionality is as of now only compatible with `dt_saving` not specified.
- `cb` correspond to callbacks function. Callbacks can be used 
to extract properties of the world at each `dt_saving` time steps of your simulation.

## Constructing the callbacks
A callback has to be of the form

```julia
cb = (names = String[], agg = Function[])
```

It is a tuple, with first value corresponding to the names of the aggregate properties of the world.
The second correspond to the aggregation functions.

We provide here an example on how to extract the ``\\gamma`` diversity of a simulation biological population. 
``\\gamma``` diversity can be calculated as the variance of the trait distribution of the population.
Here is how we write the function
```julia
cb = (names = ["gamma_div"], agg = Function[w -> var((get_x(w,1)))])
```

# Example

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
p = Dict{String,Any}();@pack! p = D,mu,NMax

myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
w0 = World(myagents, myspace, p, 0.)
w1 = copy(w0)
sim = run!(w1,Gillepsie(),tend,b,d,cb=cb,dt_saving = .1)
```


## Accessing the callbacks

You can easily access the properties, using `sim` as you would for a usual `Dictionary`.

```julia
using Plots
plot(get_tspan(sim),sim["gamma_div"])
```
"""
function run!(w::World{A,S,T},
                alg::L,
                tend::Number,
                b,
                d;
                dt_saving=nothing,
                t_saving_cb=nothing,
                cb=nothing) where {A,S,T,L<:AbstractAlg}
    # argument check
    _check_timedep(b,d)
    (!isnothing(dt_saving) && !isnothing(dt_saving)) ? ArgumentError("For now, can not specify both `dt_saving` and `t_saving_cb`") : nothing
    isnothing(dt_saving) ? dt_saving =  tend + 1. : nothing
    isnothing(t_saving_cb) ? t_saving_cb =  [tend + 1.] : nothing

    # var init
    n=length(w);
    NMax = maxsize(w)
    t = .0
    i = 1;j=1;dt = 0.
    sim = Simulation(w,cb=cb)
    if L <: Gillepsie
        update_rates!(w,alg,b,d)
    end

    # start
    while t<tend
        if dt < 0
            throw("We obtained negative time step dt = $dt at event $i")
        elseif length(w) == NMax
            @info "We have reached the maximum number of individuals allowed"
            break
        elseif length(w) == 0
            @info "All individuals have died :("
            break
        end
        if  t - get_tend(sim) >= dt_saving
            # @info "saving world @ t = $(t)/ $(tend)"
            add_entry!(sim,w,cb)
        end
        if  t >= first(t_saving_cb)
            # @info "saving callback only @ t = $(t)/ $(tend)"
            add_entry_cb_only!(sim,w,cb)
            popfirst!(t_saving_cb)
        end
        dt = updateWorld!(w,alg,b,d)
        t +=  dt
        i += 1
    end
    # Saving last time step
    add_entry!(sim,w,cb)
    @info "simulation stopped at t=$(t), after $(i) steps"
    return sim
end

"""
    _check_timedep(b,d)

Checks number of arguments of functions,
and throws error if problem
"""
function _check_timedep(b,d)
    if numargs(b) < 2
        throw(ArgumentError("Birth function needs `X` and `t` arguments"))
    end
    if numargs(d) < 3
        throw(ArgumentError("Death function needs `X`, `Y` and `t` arguments"))
    end
end
