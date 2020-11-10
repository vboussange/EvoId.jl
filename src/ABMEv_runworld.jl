"""
    function run!(w::World{A,S,T},alg::L,tend::Number,b,d;dt_saving=nothing,cb=(names = String[],agg =nothing))

Run `w` with algorithm `alg`, until `tend` is reached. User needs to provide `b` the birth function,
which takes as arguments `X,t`, and provide `d` the death function, with arguments `X,Y,t`.
Returns a `Simulation` type.
- `worldall` stores the world every `p["dt_saving"]` time steps.
If `dt_saving` not specified, `sim` contains an array of two elements,
first corresponding to initial conditions and last corresponding to world in the last time step.
- `cb` correspond to callbacks function. Look at the documentation for more information
- the run stops if the number of agents reaches`p["NMax"]`.
"""
function run!(w::World{A,S,T},alg::L,tend::Number,b,d;
                dt_saving=nothing,
                cb=(names = String[],agg =nothing)) where {A,S,T,L<:AbstractAlg}
    _check_timedep(b,d)
    n=size(w);
    NMax = maxsize(w)
    t = .0
    i = 1;j=1;dt = 0.
    isnothing(dt_saving) ? dt_saving =  tend + 1. : nothing
    sim = Simulation(w,cb=cb)
    if A <: AbstractAgent{AA,Rates{true}} where {AA}
        update_rates!(w,alg,b,d)
    end
    while t<tend
        if dt < 0
            throw("We obtained negative time step dt = $dt at event $i")
        elseif size(w) == NMax
            @info "All individuals have died :("
            break
        elseif size(w) == 0
            @info "We have reached the maximum number of individuals allowed"
            break
        end
        if  t - get_tend(sim) >= dt_saving
            @info "saving world @ t = $(t)/ $(tend)"
            add_entry!(sim,w)
        end
        dt = updateWorld!(w,alg,b,d)
        t +=  dt
        i += 1
    end
    # Saving last time step
    add_entry!(sim,w)
    @info "simulation stopped at t=$(t), after $(i) generations"
    return sim
end

"""
    function _correct_timedep!(p::Dict)

checks time dependency of birth and death functions,
and overloads the function if not provided
"""
function _check_timedep(b,d)
    if numargs(b) < 2
        throw(ArgumentError("Birth function needs `X` and `t` arguments"))
    end
    if numargs(d) < 3
        throw(ArgumentError("Death function needs `X`, `Y` and `t` arguments"))
    end
end
