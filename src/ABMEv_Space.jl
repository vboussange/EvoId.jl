
"""
    function increment_x!(a::Agent{StdAgent,U},t::U,p::Dict) where U
This function increments agent by random numbers specified in p
ONLY FOR CONTINUOUS DOMAINS
"""
function increment_x!(a::Agent{StdAgent,U},t,p::Dict) where U
    tdim = length(p["D"])
    reflected = haskey(p,"reflected") ? p["reflected"] : false
    if reflected
        inc = [get_inc_reflected(get_x(a,1),p["D"][1] *randn())]
        if  tdim > 1
            inc = vcat(inc,(rand(tdim-1) < p["mu"][2:end]) .* p["D"][2:end] .* randn(tdim-1))
        end
    else
        # inc = yes no mutation * mutation
        inc = (rand(tdim) < vec(p["mu"])) .* vec(p["D"][:]) .* randn(tdim)
    end
    a.x_history = hcat(a.x_history, get_x(a) + reshape(inc,:,1));
    push!(a.t_history,t)
 end

 """
     function increment_x!(a::Agent{MixedAgent,U},t::U,p::Dict) where U
 This function increments first trait of agent with integer values, that are automatically reflected between 1 and p["nodes"].
Other traits are incremented as usual.
TODO : make it work for a graph type landscape, where domain is not a line anymore.
 """
 function increment_x!(a::Agent{MixedAgent,U},t,p::Dict) where U
     tdim = length(p["D"])
     inc = [round(get_inc_reflected(get_x(a,1),p["D"][1] *randn(),1,p["nodes"]))]
     if  tdim > 1
         inc = vcat(inc,(rand(tdim-1) < p["mu"][2:end]) .* p["D"][2:end] .* randn(tdim-1))
     end
     a.x_history = hcat(a.x_history, get_x(a) + reshape(inc,:,1));
     push!(a.t_history,t)
end

"""
function get_inc_reflected(x::Number,inc::Number,s=-1,e=1)
    Here we increment the trajectory of trait 1 such that it follows a reflected brownian motion (1D)
"""
function get_inc_reflected(x::Number,inc::Number,s=-1,e=1)
    if x + inc < s
        inc = 2 * ( s - x ) - inc
    elseif  x + inc > e
        inc = 2 * ( e - x ) - inc
    else
        return inc
    end
    get_inc_reflected(x,inc,s,e)
end
