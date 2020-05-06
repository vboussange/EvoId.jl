using RecipesBase
using Colors
import KernelDensity:kde,pdf
@recipe function plot(world::Array{U},p;what=["x","H"],trait = 1) where U <: Union{Missing,Agent}
    tot_dim = size(world,2)*size(world,1)
    # We reduce time interval if it is too big
    if tot_dim > 1e6 && size(world,2) >= 200
        p = copy(p)
        idx_reduced = floor.(Int,range(1,size(world,2),length = 200))
        p["tspan" ] = p["tspan"][idx_reduced]
        world = world[:,idx_reduced]
    end
    if count(ismissing,world) > 0
        tspan_ar = vcat([p["tspan"][i]*ones(Int(p["NMax"] - count(ismissing,world[:,i]))) for i in 1:length(p["tspan"]) ]...);
    else
        tspan_ar = repeat(p["tspan"],inner = size(world,1))
    end
    # tspan = Float64.(tspan)
    world_sm = collect(skipmissing(world))
    if "x" in what
        d_i = []
        for i in 1:size(world,2)
            x = get_x.(skipmissing(world[:,i]),trait)
            append!(d_i,pdf(kde(x),x))
        end
        @series begin
        xarray = get_xarray(world_sm,trait)
            seriestype := :scatter
            markercolor := eth_grad_small[d_i ./ maximum(d_i)]
            # markercolor := :blue
            markerstrokewidth := 0
            alpha :=1.
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
        d_i = []
        world_df_g = groupby(world2df(world_sm),:x1)
        for world_df in world_df_g
            x = world_df.x2
            append!(d_i,pdf(kde(x),x))
        end
        # TODO: we stopped here
        @series begin
        xarray = get_xarray(world_sm,trait)
            seriestype := :scatter
            markercolor := eth_grad_small[d_i ./ maximum(d_i)]
            # markercolor := :blue
            markerstrokewidth := 0
            alpha :=1.
            xlabel := "time"
            ylabel := "trait value"
            label := ""
            grid := false
            # markersize := 2.3/1000*size(world_sm,1)
            tspan_ar[:],xarray[:]
        end
    end
    if "geo" in what
        @series begin
        xarray = get_geo.(world_sm)
            seriestype := :scatter
            markercolor := "blue"
            markerstrokewidth := 0
            alpha :=.1
            xlabel := "time"
            ylabel := "trait value"
            label := ""
            # markersize := 2.3/1000*size(world_sm,1)
            tspan_ar[:],xarray[:]
        end
    end
    if "3dgeo" in what
        @series begin
        xarray = get_geo.(world_sm)
        yarray = get_xarray(world_sm,2)
            seriestype := :scatter3d
            markercolor := "blue"
            markerstrokewidth := 0
            alpha :=.1
            xlabel := "time"
            ylabel := "geotrait"
            zlabel := "trait value"
            label := ""
            # markersize := 2.3/1000*size(world_sm,1)
            tspan_ar,xarray[:],yarray[:]
        end
    end
    if "3d" in what
        @series begin
        xarray = get_xarray(world_sm,1)
        yarray = get_xarray(world_sm,2)
            seriestype := :scatter3d
            markercolor := "blue"
            markerstrokewidth := 0
            alpha :=.1
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
