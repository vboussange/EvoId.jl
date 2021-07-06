using EVOID,Plots,UnPack


nodes = 9
dim_neutr = 30
geospace = DiscreteSegment(Int8(1),Int8(nodes))
adaptivespace = RealSpace{1,Float16}()
myspace = (geospace,adaptivespace)
sigma_K = 1.;
sigma_a = .8;
K0 = 1000;
mu = [1.,1.]
a = 1
b(X,t) = gaussian(X[2],X[1] * a,sigma_K) / nodes
d(X,Y,t) =  (X[1] ≈ Y[1]) * gaussian(X[2],Y[2],sigma_a) / K0
NMax = 2000
D = [5e-1,5e-2]
# tend = 1.5
tend = 1500
p = Dict{String,Any}();@pack! p = NMax,mu,D

myagents = [Agent(myspace,(Int8(5),Float16(5) + Float16(5e-2) * randn(Float16),),ancestors=true,rates=true) for i in 1:round(K0/nodes)]

w0 = World(myagents,myspace,p)

s = run!(w0,Gillepsie(),tend,b,d,dt_saving=5);
Plots.plot(s, ylabel = "Adaptive trait",trait=2)

using Printf
anim = @animate for i in 1:get_size(s)
    Plots.plot(s, ylabel = "Adaptive trait",
        trait=(1,2),
        time=i,
        xaxis = "Position",
        yaxis = "Adaptive trait",
        title = @sprintf("               t = %1.2f",s.tspan[i]),
        ylims = (0,11),
        xlims = (1,11))
end
gif(anim,joinpath(@__DIR__, "gradient_2Dtrait.gif"),fps = 13)

world = get_world(s,get_size(s))
shistall = get_xhist.(world[:],2)
thist = get_thist.(world[:])
aplot = Plots.plot(thist,shistall,
            linecolor = eth_grad_std[0.],
            label = "",
            # title = latexstring("\\sigma_\\mu=",@sprintf("%1.2f",world.p["D"][2]),", \\sigma_D=",@sprintf("%1.2f",world.p["D"][1])),
            grid = false,
            xlabel = "time",
            ylabel = "Lineages adaptive trait",
            ylims = (-0.5,10.5)
            )
