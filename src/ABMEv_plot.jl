using RecipesBase
using Colors
import KernelDensity:kde,pdf
"""
    function plot(world::Array{U},p;what=["x","H"],trait = 1,tplot = false) where U <: Union{Missing,Agent}

# ARGS
- `what = ["x","H"]`: the plots you want to obtain
- `trait = 1`: the trait that will plotted regarding what you asked. `trait = 0` will plot the geotrait
- `tplot = 0` used when calling xs, as it plots a snapshot of the world at a particular time
It should correspond to an integer, as it indexes the column to plot

# Options available
- `"x"`
- `"xs"`
"""

@recipe function plot(world::Array{U},p;what=["x","H"],trait = 1,tplot = 0) where U <: Union{Missing,Agent}
    tot_dim = size(world,2)*size(world,1)
    # We reduce time interval if it is too big
    # if tot_dim > 1e6 && size(world,2) >= 200
    #     p = copy(p)
    #     idx_reduced = floor.(Int,range(1,size(world,2),length = 200))
    #     p["tspan" ] = p["tspan"][idx_reduced]
    #     world = world[:,idx_reduced]
    # end
    # second condition is to make sure that the world corresponds to all the time steps
    # If not, then we can not get "x"
    if count(ismissing,world) > 0 && length(p["tspan"]) == size(world,2)
        tspan_ar = vcat([p["tspan"][i]*ones(Int(p["NMax"] - count(ismissing,world[:,i]))) for i in 1:length(p["tspan"]) ]...);
    else
        tspan_ar = repeat(p["tspan"],inner = size(world,1))
    end
    # tspan = Float64.(tspan)
    tend = p["tspan"][tplot > 0 ? tplot : length(p["tspan"])]
    world_sm = clean_world(world)
    if "x" in what
        d_i = []
        for i in 1:size(world,2)
            x = get_x(clean_world(world[:,i]),p["tspan"][i],trait)[:]
            append!(d_i,pdf(kde(x),x))
        end
        @series begin
            if length(world_sm) !== length(tspan_ar[:])
                throw(DimensionMismatch("You want to plot a world with missing data"))
            end
        xarray = get_x.(world_sm,tspan_ar[:],trait)
            seriestype := :scatter
            markercolor := eth_grad_small[d_i ./ maximum(d_i)]
            # markercolor := :blue
            markerstrokewidth := 0
            seriesalpha :=1.
            xlabel := "time"
            ylabel := "trait value"
            label := ""
            grid := false
            # markersize := 2.3/1000*size(world_sm,1)
            tspan_ar[:],xarray[:]
        end
    end
    # we use this for discrete agents
    # world should be a one dimensional vector, corresponding to one time step only
    if "xs" in what
        d_i = []; xt_array = []; x1_array = []
        world_df_all = world2df(clean_world(world[:, tplot > 0 ? tplot : size(world,2) ]),tend,true)
        world_df_g = groupby(world_df_all,:x1)
        for world_df in world_df_g
            if trait == 0
                x = Float64.(world_df.g)
            else
                # fitness occupies first spot
                x = world_df[:,trait+1] ;
            end
            x1 =  world_df.x1;
            append!(d_i,pdf(kde(x),x))
            append!(xt_array,x)
            append!(x1_array,x1)
        end
        # TODO: we stopped here
        @series begin
            seriestype := :scatter
            markercolor := eth_grad_small[d_i ./ maximum(d_i)]
            # markercolor := :blue
            markerstrokewidth := 0
            # seriesalpha := 1.
            xaxis := "geographical position"
            xticks :=  sort!(unique(world_df_all.x1))
            yaxis := "trait value"
            label := ""
            grid := false
            marker := (:rect,20,1.)
            x1_array[:],xt_array[:]
        end
    end
    if "gs" in what
        _world = clean_world(world[:, tplot > 0 ? tplot : length(p["tspan"]) ])
        y = get_x(_world,tend,2)[:]
        x = get_x(_world,tend,0)[:]
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
    if "3dgeo" in what
        d_i = []
        for i in 1:size(world,2)
            _world = clean_world(world[:,i])
            x = get_x(_world,p["tspan"][i],2)[:]
            y = get_x(_world,p["tspan"][i],0)[:]
            X = hcat(x,y)
            # d = kde(X)
            # di_temp = diag(pdf(d,X[:,1],X[:,2]))
            di_temp = y
            di_temp = (di_temp .- minimum(di_temp)) ./ (maximum(di_temp) .- minimum(di_temp))
            # here we normalise with respect to maximum value at each time step
            append!(d_i,di_temp)
        end
        @series begin
        xarray = get_geo.(world_sm,tspan_ar)
        yarray = get_x(world_sm,2)
            seriestype := :scatter3d
            markercolor := eth_grad_std[d_i ./ 1.]
            markerstrokewidth := 0
            seriesalpha :=.1
            xlabel := "time"
            ylabel := "geotrait"
            zlabel := "trait value"
            label := ""
            # markersize := 2.3/1000*size(world_sm,1)
            tspan_ar,xarray[:],yarray[:]
        end
    end
    if "3d" in what
        d_i = []
        for i in 1:size(world,2)
            x = get_x(clean_world(world[:,i]),p["tspan"][i],1)[:]
            y = get_x(clean_world(world[:,i]),p["tspan"][i],2)[:]
            X = hcat(x,y)
            d = kde(X)
            di_temp = diag(pdf(d,X[:,1],X[:,2]))
            di_temp = (di_temp .- minimum(di_temp)) ./ (maximum(di_temp) .- minimum(di_temp))
            append!(d_i,di_temp)
        end
        @series begin
        xarray = get_x(world_sm,1)
        yarray = get_x(world_sm,2)
            seriestype := :scatter3d
            markercolor := eth_grad_small[d_i ./ maximum(d_i)]
            markerstrokewidth := 0
            seriesalpha :=.1
            xlabel := "time"
            ylabel := "position"
            zlabel := "trait value"
            label := ""
            # markersize := 2.3/1000*size(world_sm,1)
            tspan_ar,xarray[:],yarray[:]
        end
    end
    # if "H" in what
    #     @series begin
    #         x = get_x.(world_sm,trait)
    #         linewidth := 2
    #         seriestype := :line
    #         label := "Interconnectedness"
    #         tspan,N^2 / 2 .* [H_discrete(x[:,i]) for i in tspan]
    #     end
    # end
    if "var" in what
        @series begin
            linewidth := 2
            seriestype := :line
            label := "Variance"
            xlabel := "Time"
            ylabel := "Variance"
            p["tspan"],var(world_sm,trait=trait)[:]
        end
    end
    if "vargeo" in what
        @series begin
            linewidth := 2
            seriestype := :line
            label := "Variance of geotrait"
            xlabel := "Time"
            ylabel := "Variance"
            p["tspan"],i->first(covgeo(world_sm[:,Int(i)]))
        end
    end
    # if "density_t" in what
    #     @series begin
    #         linewidth := 2
    #         seriestype := :plot3d
    #         label := "Variance of geotrait"
    #         xlabel := "Time"
    #         ylabel := "Variance"
    #         p["tspan"],i->first(covgeo(world_sm[:,Int(i)]))
    #     end
    # end
end
