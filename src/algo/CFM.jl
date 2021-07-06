# CFM = Champagnat Ferriere Meleard
# This contains all methods for algorithm described in article
# DOI : 10.1016/j.tpb.2005.10.004

# In this file lie the function for Gillepsie algorithm
"""
$(TYPEDEF)
"""
struct CFM <: AbstractAlg end
export CFM

function updateBirthEvent!(w::World,::CFM,mum_idx::Int)
    # updating competition only the two columns corresponding to agent idx
    offspring = give_birth(mum_idx,w)
    addAgent!(w,offspring)
end

function updateDeathEvent!(world::World,::CFM,i::Int)
    removeAgent!(world,i)
end

"""
$(SIGNATURES)
Updating rule for CFM setting.
algorithm described in article
DOI : 10.1016/j.tpb.2005.10.004
"""
function updateWorld!(w::World{A,S,T},c::CFM,b,d) where {A,S,T}
    # total update world
    @unpack bm,dm = parameters(w)
    alive = agents(w)
    # Total rate of events
    n = size(w)
    CbarI = bm + dm*(n+1)
    dt = rand(Exponential(CbarI)) / n
    update_clock!(w,dt)
    i = rand(1:n)
    x = get_x(w[i])
    W = rand()
    if dt > 0.
        deathprob = (sum(d.(get_x.(alive),Ref(x),w.t)) .- d(x,x,w.t)) / CbarI
        birthprob = b(x,w.t) / CbarI
        if W <= deathprob
            updateDeathEvent!(w,c,i)
        elseif W <= deathprob + birthprob
            updateBirthEvent!(w,c,i)
        end
        return dt
    else
        return -1
    end
end
