struct Agent
   # history of traits for geotraits
    x_history
    # birth time of ancestors
    t_history
    # death rate
    d
    #birth rate
    b
end

mutable struct mAgent
        # history of traits for geotraits
         x_history
         # birth time of ancestors
         t_history
         # death rate
         d
         #birth rate
         b
end

agentarr = [Agent([x],.1,0.,0.) for x in randn(1000)]

magentarr = [mAgent([x],.1,0.,0.) for x in randn(1000)]

using BenchmarkTools

@btime [magentarr[i].x_history[1] = randn() for i in rand(1:length(magentarr),10000)]
# 791.487 μs (14708 allocations: 386.17 KiB)
@btime [agentarr[i].x_history[1] = randn() for i in rand(1:length(agentarr),10000)]
# 968.933 μs (34728 allocations: 1.29 MiB)