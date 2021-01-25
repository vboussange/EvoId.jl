"""
    function run!(w::World{A,S,T},alg::L,tend::Number,b,d;dt_saving=nothing,cb=(names = String[],agg =nothing))

Run `w` with algorithm `alg`, until `tend` is reached. User needs to provide `b` the birth function,
which takes as arguments `X,t`, and provide `d` the death function, with arguments `X,Y,t`.
Returns a `Simulation` type.
- if `dt_saving` specified, world is saved every time steps.
If `dt_saving` not specified, `sim` contains an array of two elements,
first corresponding to initial conditions and last corresponding to world in the last time step.
- if `t_saving_cb::Array{Float64}` specified, callbacks are computed at each steps time specified in the array.
This functionality is as of now only compatible with `dt_saving` not specified.
- `cb` correspond to callbacks function. Look at the documentation for more information
- the run stops if the number of agents reaches`p["NMax"]`.
"""
function run!(w::World{A,S,T},alg::L,tend::Number,b,d;
                dt_saving=nothing,
                t_saving_cb=nothing,
                cb=nothing) where {A,S,T,L<:AbstractAlg}
    # argument check
    _check_timedep(b,d)
    (!isnothing(dt_saving) && !isnothing(dt_saving)) ? ArgumentError("For now, can not specify both `dt_saving` and `t_saving_cb`") : nothing
    isnothing(dt_saving) ? dt_saving =  tend + 1. : nothing
    isnothing(t_saving_cb) ? t_saving_cb =  [tend + 1.] : nothing

    # var init
    n=size(w);
    NMax = maxsize(w)
    t = .0
    i = 1;j=1;dt = 0.
    sim = Simulation(w,cb=cb)
    if A <: AbstractAgent{AA,Rates{true}} where {AA}
        update_rates!(w,alg,b,d)
    end

    # start
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
