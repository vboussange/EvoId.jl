cd(@__DIR__)
using Revise,ABMEv
a = 0;
sigma_K = .9;
sigma_a = 1.251;
K0 = 1000;
# K(X) = gaussian(X[1],0.,sigma_K)
K(X) = 1 - 0.125 * sum(X.^2)
α(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
# α(X,Y) = 0.
p_default = Dict(
        "alpha" => α,
        "K" => K,
        "D" => [1e-2],
        "mu" => [.1],
        "tend" => 1000.)
na_init = K0
agents = [Agent( [1e-2]  .* randn(1) .- .5) for i in 1:K0]
@time worldall,p_default["tspan"] = runWorld_store_WF(p_default,agents,reflected=false);
# ======================================================================
using JLD2
@save "wrightfisher_test.jld2" worldall p_default
using Plots
Plots.plot(worldall,p_default,what = ["var"])

var(agents)
