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
function updateWorld!(w::World{A,S,T},c::CFM) where {A,S,T}
    # total update world
    @unpack Cbar,d,b = parameters(w)
    alive = agents(w)
    # Total rate of events
    n = size(w)
    ∑ = (n + 1)
    dt = - log(rand(T))/(Cbar * ∑)
    update_clock!(w,dt)
    i = rand(1:n)
    x = get_x(w[i])
    W = rand()
    if dt > 0.
        deathprob = (sum(d.(get_x.(alive),Ref(x))) - d(x,x)) / (Cbar*(n+1))
        birthprob = b(x) / (Cbar*(n+1))
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