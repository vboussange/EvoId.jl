# for now we consider that competition is local within an array

function  update_rates_std!(world,C,p::Dict,t::Float64)
    # Competition matrix
    traits = get_x.(world)
    # traits = get_xhist.(world)
    N = length(traits)
    # C = SharedArray{Float64}((N,N))
    # Here you should do a shared array to compute in parallel
    @sync @distributed for i in 1:(N-1)
        for j in i+1:N
            C[i,j] = α(traits[i],traits[j],p["n_alpha"],p["sigma_a"])
            C[j,i] = C[i,j]
        end
    end
    # Here we can do  it in parallel as well
    for (i,a) in enumerate(world)
        # we only update death rate
        # this could be imptrove since \alpha is symmetric, by using a symmetric matrix
        a.d = sum(C[i,:]) / p["K0"]
        # /!| not ideal to assign at every time step the birth rate that is constant
        a.b = K(traits[i][:,end],1.,p["n_K"],p["sigma_K"])
    end
end

function  update_rates_graph!(world,C,p::Dict,t::Float64)
    for e in edges(p["g"])
        # agents on e
        aidx_e = findall(a -> get_x(a,1)==[e],world)
        na = length(aidx_e)
        for i in aidx_e
            world[i].d = na^2
            world[i].b = 1
        end
    end
end

function  update_rates_2D!(world,C,p::Dict,t::Float64)
    # Competition matrix
    traits = get_x.(world)
    # traits = get_xhist.(world)
    N = length(traits)
    # C = SharedArray{Float64}((N,N))
    # Here you should do a shared array to compute in parallel
    @sync @distributed for i in 1:(N-1)
        C[i,i] = 1.
        for j in i+1:N
            # be careful, for xhist what is return is an array hence no need of []
            C[i,j] = α(traits[i],traits[j],p["n_alpha"],p["sigma_a"])
            C[j,i] = C[i,j]
        end
    end
    C[N,N] = 1.
    # Here we can do  it in parallel as well
    for (i,a) in enumerate(world)
        # we only update death rate
        # a.d = sum(C[i,:]) / K(traits[i][2:2],p["K0"],p["n_K"],p["sigma_K"], μ = p["a"]*traits[i][1])
        a.d = sum(C[i,:])
        # /!| not ideal to assign at every time step the birth rate that is constant
        a.b = KK(traits[i][1:1],
                    p["K0"],
                    p["n_K"],
                    p["sigma_K"][1:1],
                    split_merge_move(t),
                    -split_merge_move(t))/K(traits[i][2:2],
                    p["K0"],
                    p["n_K"],
                    p["sigma_K"][2:2])
    end
end

function  update_rates_grad2D!(world,C,p::Dict,t::Float64)
    # Competition matrix
    traits = get_x.(world)
    # traits = get_xhist.(world)
    N = length(traits)
    # C = SharedArray{Float64}((N,N))
    # Here you should do a shared array to compute in parallel
    @sync @distributed for i in 1:(N-1)
        C[i,i] = 1.
        for j in i+1:N
            # be careful, for xhist what is return is an array hence no need of []
            C[i,j] = α(traits[i],traits[j],p["n_alpha"],p["sigma_a"])
            C[j,i] = C[i,j]
        end
    end
    C[N,N] = 1.
    # Here we can do  it in parallel as well
    for (i,a) in enumerate(world)
        # we only update death rate
        a.d = sum(C[i,:])
        a.b = K(traits[i][2:2],p["K0"],p["n_K"],p["sigma_K"], μ = p["a"]*traits[i][1])
    end
end

function  update_rates_mountain!(world,C,p::Dict,t::Float64)
    # Competition matrix
    traits = get_x.(world)
    # traits = get_xhist.(world)
    N = length(traits)
    # C = SharedArray{Float64}((N,N))
    # Here you should do a shared array to compute in parallel
    @sync @distributed for i in 1:(N-1)
        C[i,i] = 1.
        for j in i+1:N
            # be careful, for xhist what is return is an array hence no need of []
            C[i,j] = α(traits[i],traits[j],p["n_alpha"],p["sigma_a"])
            C[j,i] = C[i,j]
        end
    end
    C[N,N] = 1.
    # Here we can do  it in parallel as well
    for (i,a) in enumerate(world)
        # we only update death rate
        a.d = sum(C[i,:])
        a.b = KK(traits[i][1:1],
                p["K0"],
                p["n_K"],
                p["sigma_K"][1:1],
                split_move(t),
                -split_move(t))/K(traits[i][2:2],
                p["K0"],
                p["n_K"],
                p["sigma_K"][2:2],
                μ=(tin(traits[i][1],-1.,0.)*(traits[i][1]+1) - tin(traits[i][1],0.,1.)*(traits[i][1]-1) ) * split_move(t))
    end
end

function  update_rates_std_split!(world,C,p::Dict,t::Float64)
    # Competition matrix
    traits = get_x.(world)
    # traits = get_xhist.(world)
    N = length(traits)
    # C = SharedArray{Float64}((N,N))
    # Here you should do a shared array to compute in parallel
    @sync @distributed for i in 1:(N-1)
        C[i,i] = 1.
        for j in i+1:N
            C[i,j] = α(traits[i],traits[j],p["n_alpha"],p["sigma_a"])
            C[j,i] = C[i,j]
        end
    end
    C[N,N] = 1.
    # Here we can do  it in parallel as well
    for (i,a) in enumerate(world)
        # we only update death rate
        # this could be imptrove since \alpha is symmetric, by using a symmetric matrix
        # a.d = sum(C[i,:]) / KK(traits[i][:,end],p["K0"],p["n_K"],p["sigma_K"],split_merge_move(t),-split_merge_move(t))
        a.d = sum(C[i,:])
        # /!| not ideal to assign at every time step the birth rate that is constant
        a.b =  KK(traits[i][:,end],p["K0"],p["n_K"],p["sigma_K"],split_move(t),-split_move(t))
    end
end
"""
    function runWorld_store_WF(p,world0;init = ([.0],),reflected=false,mode="std")
Wright Fisher process. Returns an array worldall with every agents.
"""
function runWorld_store_WF(p,world0;init = ([.0],),reflected=false,mode="std")
    worldall = repeat(world0,inner = (1,length(1:Int(p["tend"]))))
    N=length(world0);
    C = SharedArray{Float64}((N,N));
    newworld = copy.(world0)
    if mode == "std"
        update_rates! = update_rates_std!
    elseif mode == "2D"
        update_rates! = update_rates_2D!
    elseif mode == "grad2D"
        update_rates! = update_rates_grad2D!
    elseif mode == "mountain"
        update_rates! = update_rates_mountain!
    elseif mode == "split"
        update_rates! = update_rates_std_split!
    elseif mode == "graph"
        update_rates! = update_rates_graph!
    else
        error("Mode $mode is not recognized")
    end
    for i in 1:(Int(p["tend"])-1)
        # we have to take views, otherwise it does not affect worldall
        world = @view worldall[:,i];newworld = @view worldall[:,i+1];
        updateWorld_WF!(world,newworld ,C,p,update_rates!,Float64(i),reflected)
    end
    return worldall
end
"""
    function runWorld_store_G(p,world0;init = ([.0],),reflected=false)
Gillepsie process. Returns a tuple worldall,tspanarray
"""
function runWorld_store_G(p,world0;init = ([.0],),reflected=false)
    # we store the value of world every 100 changes by default
    tspan = zeros(1)
    i = 1;j=1;
    N=length(world0);
    tspanarray = [];
    ninit = Int(length(world0) - count(ismissing,world0))
    # Array that stores the agents.
    # We decide that we store agents ninit events where ninit is the initial population
    # length of worldall should be changed at some point
    worldall = reshape(copy.(world0),N,1)
    # we instantiate C as the biggest size it can take
    C = SharedArray{Float64}((N,N))
    update_rates_std!(skipmissing(world0),C,p,0.)
    while tspan[i]<p["tend"] &&  count(ismissing,world0) < p["NMax"] && count(ismissing,world0) > 0
        # we save every ninit times
        if mod(i,ninit) == 1
            @info "saving world @ t = $(tspan[i])/ $(p["tend"])"
            j+=1;sw = size(worldall,2);
            if sw < j
                # we double the size of worldall
                worldall = hcat(worldall,Array{Missing}(missing,N,sw))
            end
            worldall[1:Int(N - count(ismissing,world0)),j] .= copy.(collect(skipmissing(world0)));
            push!(tspanarray,tspan[i])
        end
        push!(tspan, tspan[end] + updateWorld_G!(world0,C,p,update_rates_std!,tspan,reflected))
        i += 1
    end
    return worldall[:,1:j],tspanarray
end
