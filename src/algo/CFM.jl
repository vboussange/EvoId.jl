# CFM = Champagnat Ferriere Meleard
# This contains all methods for algorithm described in article
# DOI : 10.1016/j.tpb.2005.10.004

# In this file lie the function for Gillepsie algorithm
"""
$(TYPEDEF)

Champagnat Ferriere Méléard algorithm described in 
[Champagnat and Ferriere founding article](https://linkinghub.elsevier.com/retrieve/pii/S0040580905001632). 
This algorithm is similar to `Gillepsie`, excepts that it runs faster for higher number of individuals.

Indeed, at every time step, only the fitness of the individual picked at random is evaluated. 
In order to use it, you need to feed to the dictionary parameters 
`p` a constant `Cbar<:Real` that is 
the upperbound of the maximum of the sum of the birth and death rates (cf article).

# Example
```julia
using EvoId,UnPack,Plots
myspace = (RealSpace{1,Float64}(),)
σ_b = .9;
σ_d = .7;
K0 = 1000
b(X,t) = 1.
d(X,Y,t) = gaussian(X[1],Y[1],σ_d)/K0 / gaussian(X[1],0.,σ_b)
Cbar = b([0],0.) + d([0],[0],0.)
D = (1e-2,)
mu = [.1]
NMax = 2000
tend = 1500
p = Dict{String,Any}();@pack! p = D,mu,NMax,Cbar
myagents = [Agent(myspace,(1e-2 * randn(),)) for i in 1:K0]
w0 = World(myagents,myspace,p,0.)
@time sim = run!(w0,CFM(),tend, b, d, dt_saving = 4)
```

!!! warning "Development"
    CFM gives an approximate time step. As of now, we do not manage to obtain qualitatively the same results as the Gillepsie algorithm.


Algorithm described in article
DOI : 10.1016/j.tpb.2005.10.004
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
