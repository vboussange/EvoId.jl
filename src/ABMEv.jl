__precompile__()

module ABMEv
    using Distributions,LinearAlgebra,Reexport,StatsBase
    using LightGraphs

    include("ABMEv_Agent.jl")
    include("ABMEv_WF.jl")
    include("ABMEv_Gillepsie.jl")
    include("ABMEv_runworld.jl")
    include("ABMEv_metrics.jl")
    include("ABMEv_plot.jl")
    include("ABMEv_utils.jl")


    @reexport using Distributions, DataFrames
    export update_rates!
    export MixedAgent,StdAgent,Agent,get_fitness,get_x,get_dim,get_nancestors,get_xarray,get_xhist,
        get_geo,get_b,get_d,increment_x!,get_inc_reflected,world2df,
        split_move,split_merge_move,tin,new_world_G
    export copy,runWorld_store_WF,runWorld_store_G,clean_world #,runWorld_G!,runWorld_WF!,
    export H_discrete,findclusters,var,covgeo,hamming,get_beta_div, get_alpha_div
    export update_afterbirth_std!,update_afterdeath_std!
    export generalised_gaussian,gaussian,ma,eth_grad_std
end
