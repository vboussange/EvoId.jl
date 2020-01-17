using Plots
# function designed for Wright Fisher process
@recipe function plot(world::Array{Agent{T}},p;what=["x","H"],trait = 1) where T
    # here  we can take xhist because every agent is updated at the same time
    # world should correspond to a one dimensional array
    N = size(world,1)
    tspan = collect(1:Int(p["tend"]))
    if "x" in what
        @series begin
        xarray = get_xarray(world,trait)
            seriestype := :scatter
            markercolor := "blue"
            markerstrokewidth := 0
            alpha :=.1
            xlabel := "time"
            ylabel := "trait value"
            label := ""
            markersize := 2.3/1000*size(world,1)
            repeat(tspan,inner = N),xarray[:]
        end
    end
    if "geo" in what
        @series begin
        xarray = get_geo.(world)
            seriestype := :scatter
            markercolor := "blue"
            markerstrokewidth := 0
            alpha :=.1
            xlabel := "time"
            ylabel := "trait value"
            label := ""
            markersize := 2.3/1000*size(world,1)
            repeat(tspan,inner = N),xarray[:]
        end
    end
    if "3dgeo" in what
        @series begin
        xarray = get_geo.(world)
        yarray = get_xarray(world,2)
            seriestype := :scatter3d
            markercolor := "blue"
            markerstrokewidth := 0
            alpha :=.1
            xlabel := "time"
            ylabel := "geotrait"
            zlabel := "trait value"
            label := ""
            markersize := 2.3/1000*size(world,1)
            repeat(tspan,inner = N),xarray[:],yarray[:]
        end
    end
    if "3d" in what
        @series begin
        xarray = get_xarray(world,1)
        yarray = get_xarray(world,2)
            seriestype := :scatter3d
            markercolor := "blue"
            markerstrokewidth := 0
            alpha :=.1
            xlabel := "time"
            ylabel := "position"
            zlabel := "trait value"
            label := ""
            markersize := 2.3/1000*size(world,1)
            repeat(tspan,inner = N),xarray[:],yarray[:]
        end
    end
    # if "H" in what
    #     @series begin
    #         x = get_x.(world,trait)
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
            tspan,var(world,trait=trait)[:]
        end
    end
    if "vargeo" in what
        @series begin
            linewidth := 2
            seriestype := :line
            label := "Variance of geotrait"
            xlabel := "Time"
            ylabel := "Variance"
            tspan,i->first(covgeo(world[:,Int(i)]))
        end
    end
end

@recipe function plot(world::Array{Union{Agent{T},Missing}},p;what=["x","H"],trait = 1) where T
    # here world should correspond to a multidimensional array
    xarray = vcat(get_x.(skipmissing(world),trait)...)
    N = size(xarray,1)
    tspan = p["tspan"]
    tspanarray = vcat([tspan[i]*ones(Int(p["NMax"] - count(ismissing,world[:,i]))) for i in 1:length(tspan) ]...)
    if "x" in what
        @series begin
            seriestype := :scatter
            markercolor := "blue"
            markerstrokewidth := 0
            alpha :=.1
            xlabel := "time"
            ylabel := "trait value"
            label := ""
            markersize := 2.3/1000*size(world,1)
            tspanarray,xarray
        end
    end
    if "H" in what
        @series begin
            x = get_x.(world,trait)
            linewidth := 2
            seriestype := :line
            label := "Interconnectedness"
            tspan,N^2 / 2 .* [H_discrete(x[:,i]) for i in tspan]
        end
    end
end
