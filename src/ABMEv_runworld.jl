function  update_rates_graph!(world,C,p::Dict,t::Float64)
    for e in edges(p["g"])
        # agents on e
        aidx_e = findall(a -> get_x(a,1)==e,world)
        na = length(aidx_e)
        for i in aidx_e
            world[i].d = na^2
            world[i].b = 1
        end
    end
end

"""
$(SIGNATURES)
Run `w` with algorithm `alg`, until `tend` is reached.
Returns a `Simulation` type.
- `worldall` stores the world every `p["dt_saving"]` time steps.
If `dt_saving` not specified, `sim` contains an array of two elements,
first corresponding to initial conditions and last corresponding to world in the last time step.
>:warning: if you choose `nagents = 1` then nothing is saved until the end of the simulation.
"""
# function run(w::World{AbstractAgent{A,R},S,T},g::G;dt_saving=nothing,callbacks=nothing) where {G<:Gillepsie,A,R,S,T}
function run!(w::World{A,S,T},alg::L,tend::Number;
                dt_saving=nothing,
                cb=(names = String[],agg =nothing)) where {A,S,T,L<:AbstractAlg}
    n=size(w);
    NMax = maxsize(w)
    t = .0
    i = 1;j=1;dt = 0.
    isnothing(dt_saving) ? dt_saving =  tend + 1. : nothing
    sim = Simulation(w,cb=cb)
    if A <: AbstractAgent{AA,Rates{true}} where {AA}
        update_rates!(w,alg)
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
        dt = updateWorld!(w,alg)
        t +=  dt
        i += 1
    end
    # Saving last time step
    add_entry!(sim,w)
    @info "simulation stopped at t=$(t), after $(i) generations"
    return sim
end
