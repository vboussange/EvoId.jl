__precompile__()

module ABMEv
    using Distributions,LinearAlgebra,Reexport,StatsBase
    using LightGraphs
    using UnPack
    using DocStringExtensions

    include("ABMEv_Space.jl")
    include("ABMEv_Agent.jl")
    include("ABMEv_world.jl")
    include("ABMEv_Sim.jl")
    include("ABMEv_metrics.jl")
    # include("ABMEv_plot.jl")
    include("ABMEv_utils.jl")
    include("algo/ABMEv_WF.jl")
    include("algo/ABMEv_Gillepsie.jl")
    include("algo/ABMEv_CFM.jl")
    include("ABMEv_runworld.jl")


    @reexport using Distributions, DataFrames

    export GraphSpace,ContinuousSegment,DiscreteSegment,RealSpace,
        AbstractSpacesTuple,get_inc
    export update_rates!
    export AbstractAgent,Agent,get_fitness,get_x,get_t,get_dim,
        nancestors,get_xarray,get_xhist,
        get_thist,get_geo,get_b,get_d,increment_x!,get_inc_reflected,world2df,
        split_move,split_merge_move,tin,new_world_G
    export World,parameters,time,space,agents,size,maxsize,addAgent!,removeAgent!
    export run!,give_birth,updateWorld!,update_clock!,updateBirthEvent!,
        updateDeathEvent!#,runWorld_G!,runWorld_WF!,
    export Simulation,add_entry!,get_tend,get_size,get_tspan
    export H_discrete,findclusters,var,covgeo,hamming,get_beta_div, get_alpha_div,
        get_dist_hist,get_pairwise_average_isolation,
        get_local_pairwise_average_isolation,
        truncvar,get_xhist_mat
    export update_afterbirth_std!,update_afterdeath_std!
    export generalised_gaussian,gaussian,ma,geomsmooth,arithsmooth,eth_grad_std,
        DiversityFunction,geomsmooth2D,arithsmooth2D,interpolate_df,groupby
end
