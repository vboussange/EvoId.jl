using ABMEv,LightGraphs,UnPack

##### Genotype space#####
nodes = 9
dim_neutr = 1000
magicprop = 523728 / 32896
g = SimpleGraph{Int16}(dim_neutr,round(Int16,dim_neutr * magicprop))
initnode = argmax(eigenvector_centrality(g)) # This is the central node the we will use to instantiate the populations
myspace = (DiscreteSegment(Int8(1),Int8(nodes)),GraphSpace(g)) # union of vector spaces

K0 = 1000
mu = [1.,1.]
b(X,t) = 1 / nodes
d(X,Y,t) = (X[1] ≈ Y[1]) / K0
D = (3e-1,5e-1)

NMax = 2000
# tend = 1.5
tend = 200
p_default = Dict{String,Any}();@pack! p_default = NMax,mu,D
myagents = [Agent(myspace,(Int8(5),initnode),ancestors=true,rates=true) for i in 1:round(K0/nodes)]
w0 = World(myagents,myspace,p_default)

@time sim = run!(w0,Gillepsie(),tend,b,d,dt_saving=3.)

using GraphPlot,StatsBase

## Plotting genetic space
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
d_i = [ (_d_i .- minimum(_d_i)) ./ (maximum(_d_i) .- minimum(_d_i)) for _d_i in d_i ]
using Printf,Cairo,Compose
for i in 1:length(d_i)
        gp =gplot(g,locs_x,locs_y,
                nodefillc=eth_grad_std[d_i][i])
        draw(PNG(joinpath(@__DIR__,"gen_struct/gif_genetic_structure_$i.png"), 16cm, 16cm), gp)

        # Plots.plot(xarray,p_default,what=["gs"],trait=2,tplot=t,
        # title = @sprintf("               t = %1.2f",p_default["tspan"][t]),
        # ylims = (1,p_default["nodes"]),
        # xlims = (0,maximum(xgeo))
        # )
end
## Plotting spatial space
using Plots
# This is to plot with consistency
xarray,tspan = get_xnt(sim, trait=1)
# here we compute the number of individuals per genome, throught all the time steps
d_i = []
for x in xarray
        _d_i = zeros(nv(g))
        c = countmap(x)
        _d_i[collect(keys(c))] .= values(c)
        push!(d_i,_d_i)
end
# Here we normalize with respect to the max size of the whole world through time
d_i = [ (_d_i .- minimum(_d_i)) ./ (maximum(_d_i) .- minimum(_d_i)) for _d_i in d_i ]
anim  = @animate for i in 1:length(d_i)
        Plots.scatter(collect(1:nodes),ones(nodes),
                markercolor=ABMEv.eth_grad_small[d_i][i],
                markersize = 40,
                grid = false,
                xaxis = false,
                yaxis = false,
                label = "",
                markerstrokewidth = 0.,
                ylims = (-8,10),
                xlims = (0,10)
                )
        # Plots.plot(xarray,p_default,what=["gs"],trait=2,tplot=t,
        # title = @sprintf("               t = %1.2f",p_default["tspan"][t]),
        # ylims = (1,p_default["nodes"]),
        # xlims = (0,maximum(xgeo))
        # )
end
gif(anim,joinpath(@__DIR__, "space_genetic_struct.gif"),fps = 13)
