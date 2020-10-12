using RecipesBase
using Colors
import KernelDensity:kde,pdf
"""
    function plot(sim::Simulation;trait = 1)
Plot recipe for ABMEv.jl
# ARGS
- if `length(trait) == 1` then we scatter plot `trait` along time
- if `2 <= length(trait) <= 3` then we project world of the
last time step in the two  (three) dimensional trait space define by `trait`

!!! warning "To be implemented"
    We might want to get a 3 dimensional scatter plot
    with time, trait1 and trait2 as axis
"""

@recipe function plot(sim::Simulation;trait = 1,time = nothing)
    # world = sim.agentarray
    # tspan = sim.tspan
    # p = sim.p
    if length(trait) == 1
        xarray,tspan = get_xnt(sim, trait=trait)
        d_i = []
        for x in xarray
            push!(d_i,pdf(kde(x),x))
        end
        # Here we normalize with respect to the max size of the whole world through time
        d_i = [ (_d_i .- minimum(_d_i)) ./ (maximum(_d_i) .- minimum(_d_i)) for _d_i in d_i ]
        @series begin
            seriestype := :scatter
            markercolor := eth_grad_std[vcat(d_i...)]
            markerstrokewidth := 0
            seriesalpha := .2
            # xlabel := "time"
            # ylabel := "trait value"
            label := ""
            xlabel := "time"
            grid := false
            # markersize := 2.3/1000*size(world_sm,1)
            vcat(tspan...),vcat(xarray...)
        end
    end
    if length(trait) == 2
        isnothing(time) ? w = get_world(sim[end]) :  w = get_world(sim[time])
        y = get_x(w,trait[1]);
        x=get_x(w,trait[2])
        X = hcat(x,y)
        d = kde(X)
        # by density
        d_i = diag(pdf(d,X[:,1],X[:,2]))
        # by value
        # d_i = y
        d_i = (d_i .- minimum(d_i)) ./ (maximum(d_i) .- minimum(d_i))
        # TODO: we stopped here
        @series begin
            seriestype := :scatter
            markercolor := eth_grad_small[d_i]
            # markercolor := :blue
            markerstrokewidth := 0
            # seriesalpha := 1.
            xaxis := "geotrait"
            yaxis := "trait value"
            label := ""
            grid := false
            # marker := (:rect,20,1.)
            x,y
        end
    end
    if length(trait) == 2
        throw(ArgumentError("Plot for three traits not yet implemented"))
    end
end


import Plots:cgrad
# asymmetry towards red, blue is only a fifth of the color range
const eth_grad_small = cgrad([colorant"#1F407A", RGB(0.671,0.851,0.914),RGB(1.0,1.0,0.749), RGB(0.992,0.682,0.38),RGB(0.647,0.0,0.149),],[.0,.1])
# symmetry between red and blue
const eth_grad_std = cgrad([colorant"#1F407A", RGB(0.671,0.851,0.914),RGB(1.0,1.0,0.749), RGB(0.992,0.682,0.38),RGB(0.647,0.0,0.149),],[.0,1.])


########## OLD PLOT RECIPES ##############
# To get some inspiration

# we use this for discrete agents
# world should be a one dimensional vector, corresponding to one time step only
# if "xs" in what
#     d_i = []; xt_array = []; x1_array = []
#     world_df_all = world2df(clean_world(world[:, tplot > 0 ? tplot : size(world,2) ]),tend,true)
#     world_df_g = groupby(world_df_all,:x1)
#     for world_df in world_df_g
#         if trait == 0
#             x = Float64.(world_df.g)
#         else
#             # fitness occupies first spot
#             x = world_df[:,trait+1] ;
#         end
#         x1 =  world_df.x1;
#         append!(d_i,pdf(kde(x),x))
#         append!(xt_array,x)
#         append!(x1_array,x1)
#     end
#     @series begin
#         seriestype := :scatter
#         markercolor := eth_grad_small[d_i ./ maximum(d_i)]
#         # markercolor := :blue
#         markerstrokewidth := 0
#         # seriesalpha := 1.
#         xaxis := "geographical position"
#         xticks :=  sort!(unique(world_df_all.x1))
#         yaxis := "trait value"
#         label := ""
#         grid := false
#         # marker := (:rect,20,1.)
#         x1_array[:],xt_array[:]
#     end
# end
##

# if "3dgeo" in what
#     d_i = []
#     for i in 1:size(world,2)
#         _world = clean_world(world[:,i])
#         x = get_x(_world,tspan[i],2)[:]
#         y = get_x(_world,tspan[i],0)[:]
#         X = hcat(x,y)
#         # d = kde(X)
#         # di_temp = diag(pdf(d,X[:,1],X[:,2]))
#         di_temp = y
#         di_temp = (di_temp .- minimum(di_temp)) ./ (maximum(di_temp) .- minimum(di_temp))
#         # here we normalise with respect to maximum value at each time step
#         append!(d_i,di_temp)
#     end
#     @series begin
#     xarray = get_geo.(world_sm,tspan_ar)
#     yarray = get_x(world_sm,2)
#         seriestype := :scatter3d
#         markercolor := eth_grad_std[d_i ./ 1.]
#         markerstrokewidth := 0
#         seriesalpha :=.1
#         xlabel := "time"
#         ylabel := "geotrait"
#         zlabel := "trait value"
#         label := ""
#         # markersize := 2.3/1000*size(world_sm,1)
#         tspan_ar,xarray[:],yarray[:]
#     end
# end
# if "3d" in what
#     d_i = []
#     for i in 1:size(world,2)
#         x = get_x(clean_world(world[:,i]),tspan[i],1)[:]
#         y = get_x(clean_world(world[:,i]),tspan[i],2)[:]
#         X = hcat(x,y)
#         d = kde(X)
#         di_temp = diag(pdf(d,X[:,1],X[:,2]))
#         di_temp = (di_temp .- minimum(di_temp)) ./ (maximum(di_temp) .- minimum(di_temp))
#         append!(d_i,di_temp)
#     end
#     @series begin
#     xarray = get_x(world_sm,1)
#     yarray = get_x(world_sm,2)
#         seriestype := :scatter3d
#         markercolor := eth_grad_small[d_i ./ maximum(d_i)]
#         markerstrokewidth := 0
#         seriesalpha :=.1
#         xlabel := "time"
#         ylabel := "position"
#         zlabel := "trait value"
#         label := ""
#         # markersize := 2.3/1000*size(world_sm,1)
#         tspan_ar,xarray[:],yarray[:]
#     end
# end
# if "H" in what
#     @series begin
#         x = get_x.(world_sm,trait)
#         linewidth := 2
#         seriestype := :line
#         label := "Interconnectedness"
#         tspan,N^2 / 2 .* [H_discrete(x[:,i]) for i in tspan]
#     end
# end
# if "var" in what
#     @series begin
#         linewidth := 2
#         seriestype := :line
#         label := "Variance"
#         xlabel := "Time"
#         ylabel := "Variance"
#         tspan,var(world_sm,trait=trait)[:]
#     end
# end
# if "vargeo" in what
#     @series begin
#         linewidth := 2
#         seriestype := :line
#         label := "Variance of geotrait"
#         xlabel := "Time"
#         ylabel := "Variance"
#         tspan,i->first(covgeo(world_sm[:,Int(i)]))
#     end
# end
# if "density_t" in what
#     @series begin
#         linewidth := 2
#         seriestype := :plot3d
#         label := "Variance of geotrait"
#         xlabel := "Time"
#         ylabel := "Variance"
#         tspan,i->first(covgeo(world_sm[:,Int(i)]))
#     end
# end
