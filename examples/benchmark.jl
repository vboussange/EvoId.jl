using ABMEv,UnPack

println("--------------------------------")
println("""TESTING GILLEPSIE with popoulation""")
println("--------------------------------")

for K0 in [10,50,100,500,1000,5000]
    myspace = (RealSpace{1,Float64}(),)
    sigma_K = .9;
    sigma_a = .7;
    b(X,t) = gaussian(X[1],0.,sigma_K)
    d(X,Y,t) = gaussian(X[1],Y[1],sigma_a)/K0
    D = (1e-2,)
    mu = [.1]
    NMax = 10000
    tend = 2
    p = Dict{String,Any}();@pack! p = d,b,D,mu,NMax
    myagents = [Agent(myspace,(0,),ancestors=false,rates=true) for i in 1:K0]
    w0 = World(myagents,myspace,p,0.)
    @show K0
    @time sim = run!(w0,Gillepsie(),tend,b,d)
end

println("--------------------------------")
println("""TESTING CFM with population""")
println("--------------------------------")
for K0 in [10,50,100,500,1000,5000]
    myspace = (RealSpace{1,Float64}(),)
    sigma_K = .9;
    sigma_a = .7;
    b(X,t) = gaussian(X[1],0.,sigma_K)
    d(X,Y,t) = gaussian(X[1],Y[1],sigma_a)/K0
    D = (1e-2,)
    mu = [.1]
    NMax = 10000
    tend = 2
    dm = d([0],[0],0.);bm = 1.
    p = Dict{String,Any}();@pack! p = dm,bm,D,mu,NMax,Cbar
    myagents = [Agent(myspace,(0,),ancestors=false) for i in 1:K0]
    w0 = World(myagents,myspace,p,0.)
    @show K0
    @time sim = run!(w0,CFM(),tend,b,d)
end

println("--------------------------------")
println("""TESTING CFM with population without ancestors""")
println("--------------------------------")
for K0 in [10,50,100,500,1000,5000]
    myspace = (RealSpace{1,Float64}(),)
    sigma_K = .9;
    sigma_a = .7;
    b(X) = gaussian(X[1],0.,sigma_K)
    d(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
    D = (1e-2,)
    mu = [.1]
    NMax = 10000
    tend = 2
    Cbar=2
    p = Dict{String,Any}();@pack! p = D,mu,NMax,Cbar
    myagents = [Agent(myspace,(0,),ancestors=true) for i in 1:K0]
    w0 = World(myagents,myspace,p,0.)
    @show K0
    @time sim = run!(w0,CFM(),tend,b,d)
end

println("--------------------------------")
println("""TESTING CFM with time""")
println("--------------------------------")

for tend in [1,10,50,100]
    K0 = 1000
    myspace = (RealSpace{1,Float64}(),)
    sigma_K = .9;
    sigma_a = .7;
    b(X) = gaussian(X[1],0.,sigma_K)
    d(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
    D = (1e-2,)
    mu = [.1]
    NMax = 10000
    Cbar=2
    p = Dict{String,Any}();@pack! p = D,mu,NMax,Cbar
    myagents = [Agent(myspace,(0,),ancestors=true) for i in 1:K0]
    w0 = World(myagents,myspace,p,0.)
    @show tend
    @time sim = run!(w0,CFM(),tend,b,d)
end

println("--------------------------------")
println("""TESTING CFM with time without Ancestors""")
println("--------------------------------")

for tend in [1,10,50,100]
    K0 = 1000
    myspace = (RealSpace{1,Float64}(),)
    sigma_K = .9;
    sigma_a = .7;
    b(X) = gaussian(X[1],0.,sigma_K)
    d(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
    D = (1e-2,)
    mu = [.1]
    NMax = 10000
    Cbar=2
    p = Dict{String,Any}();@pack! p = D,mu,NMax,Cbar
    myagents = [Agent(myspace,(0,)) for i in 1:K0]
    w0 = World(myagents,myspace,p,0.)
    @show tend
    @time sim = run!(w0,CFM(),tend,b,d)
end

println("--------------------------------")
println("""TESTING Gillepsie with time""")
println("--------------------------------")
for tend in [1,10,50,100]
    K0 = 1000
    myspace = (RealSpace{1,Float64}(),)
    sigma_K = .9;
    sigma_a = .7;
    b(X) = gaussian(X[1],0.,sigma_K)
    d(X,Y) = gaussian(X[1],Y[1],sigma_a)/K0
    D = (1e-2,)
    mu = [.1]
    NMax = 10000
    # Cbar=2
    p = Dict{String,Any}();@pack! p = D,mu,NMax #,Cbar
    myagents = [Agent(myspace,(0,),ancestors=true,rates=true) for i in 1:K0]
    w0 = World(myagents,myspace,p,0.)
    @show tend
    @time sim = run!(w0,Gillepsie(),tend,b,d)
end
