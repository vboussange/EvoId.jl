using ABMEv,UnPack

##### Genotype space#####
dim_neutr = 1000
magicprop = 523728 / 32896
g = SimpleGraph{Int16}(dim_neutr,round(Int16,dim_neutr * magicprop))
initnode = argmax(eigenvector_centrality(g)) # This is the central node the we will use to instantiate the populations
myspace = (DiscreteSegment(Int8(1),Int8(nodes)),GraphSpace(g)) # union of vector spaces

K0 = 1000
mu = [1.,1.]
b(X,t) = 1 / nodes
d(X,Y,t) = (X[1] ≈ Y[1]) / K0
D = (5e-1,1.4825)

NMax = 2000
# tend = 1.5
tend = 100
p_default = Dict{String,Any}();@pack! p_default = d,b,NMax,mu,D
myagents = [Agent(myspace,(rand(Int8(1):Int8(nodes)),initnode),ancestors=true,rates=true) for i in 1:round(K0/nodes)]
w0 = World(myagents,myspace,p_default,0.)

@time sim = run!(w0,Gillepsie(),tend,dt_saving=3.)

using GraphPlot,StatsBase

# Plotting genetic space
# This is to plot with consistency
locs_x, locs_y = spring_layout(g)
xarray,tspan = get_xnt(sim, trait=2)
# here we compute the number of individuals per genome, throught all the time steps
d_i = []
for x in xarray
        _d_i = zeros(nv(g))
        c = countmap(x)
        _d_i[collect(keys(c))] .= values(c)
        push!(d_i,_d_i)
end
# Here we normalize with respect to the max size of the whole world through time
d_i = [ (_d_i .- minimum(_d_i)) ./ (maximum(_d_i) .- minimum(_d_i)) for _d_i in d_i ]anim = @animate for t in 1:size(worldall,2)
Plots.plot(worldall,p_default,what=["gs"],trait=2,tplot=t,
title = @sprintf("               t = %1.2f",p_default["tspan"][t]),
ylims = (1,p_default["nodes"]),
xlims = (0,maximum(xgeo))
)
                end
gif(anim,"gif_d2=$(p_default["D"][2])_d1=$(p_default["D"][1])_tend=$(p_default["tend"])_mu=1e0_1e0.gif",fps = 13)
