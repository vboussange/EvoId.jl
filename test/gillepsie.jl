using ABMEv
a = 0;
sigma_K = .9;
sigma_a = .7;
K(X) = gaussian(X[1],0.,sigma_K)
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/1000
# α(X,Y) = 0.
p_default = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2], # we let mutation constant and =1e-3
        "mu" => [.1], # we consider there is a probability of 1 to have mutation
        "tend" => 150.,
        "NMax" => Int(10000))
na_init = 150
world0 = new_world_G(na_init,p_default,spread = .01, offset = -.25)
@time worldall,p_default["tspan"] = runWorld_store_G(p_default,world0)
# ======================================================================
using JLD2
@save "gillepsie_test.jld2" worldall,p_default
using Plots
Plots.plot(worldall,p_default)
