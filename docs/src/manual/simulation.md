# Simulation
A `Simulation` object is returned by the function `run!`. It is a container for snapshots of the world at every `dt_saving` time steps. It renders post processing easier, through dedicated methods to obtain time series of quantities.


!!! note "Default behaviour"
    If `df_saving` is not provided, initial and last time steps will be saved.




```@autodocs
Modules = [EvoId]
Pages   = ["Sim.jl"]
```
