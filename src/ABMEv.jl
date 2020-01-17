__precompile__()

module ABMEv
    using Distributions,LinearAlgebra,Reexport,StatsBase
    using SharedArrays,Distributed,LightGraphs

    include("ABMEv_Agent.jl")
    include("ABMEv_WF.jl")
    include("ABMEv_Gillepsie.jl")
    include("ABMEv_runworld.jl")
    include("ABMEv_metrics.jl")
    include("ABMEv_plot.jl")


    @reexport using Distributions,Plots
    export update_rates!
    export Agent,get_fitness,get_x,get_xarray,get_xhist,
        get_geo,get_b,get_d,increment_x!,get_inc_reflected,
        split_move,split_merge_move,KK,tin
    export copy,runWorld_store_WF,runWorld_store_G #,runWorld_G!,runWorld_WF!,
    export H_discrete,findclusters,var,covgeo,hamming
end
