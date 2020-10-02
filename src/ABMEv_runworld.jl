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
    function runWorld_store_WF(p,world0;mode="std")
Wright Fisher process. Returns an array worldall with every agents.
"""
function runWorld_store_WF(p,world0;mode="std")
    # worldall = repeat(world0,inner = (1,length(1:Int(p["tend"]))))
    # N=length(world0);
    # newworld = copy.(world0)
    # if mode == "std"
    #     update_rates! = update_rates_std!
    # elseif mode == "2D"
    #     update_rates! = update_rates_2D!
    # elseif mode == "grad2D"
    #     update_rates! = update_rates_grad2D!
    # elseif mode == "mountain"
    #     update_rates! = update_rates_mountain!
    # elseif mode == "split"
    #     update_rates! = update_rates_std_split!
    # elseif mode == "graph"
    #     update_rates! = update_rates_graph!
    # else
    #     error("Mode $mode is not recognized")
    # end
    # for i in 1:(Int(p["tend"])-1)
    #     # we have to take views, otherwise it does not affect worldall
    #     world = @view worldall[:,i];newworld = @view worldall[:,i+1];
    #     updateWorld_WF!(world,newworld,p,update_rates!,Float64(i))
    # end
    # return worldall,collect(0:Int(p["tend"]-1))
end
"""
    function runWorld_store_G(p,world0)
Gillepsie process.
- returns a tuple 'worldall, tspanarray'
- `worldall` stores the world every `p["dt_saving"]` time steps.
If `p["dt_saving"]` not specified, it returns an array with two columns,
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
