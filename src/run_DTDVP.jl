function run_DTDVP(dt, tmax, A, H, prec; obs=[], effects=false, error=false, timed=false, savebonddims=false, Dplusmax=nothing, Dlim=50, reduceddensity=false, timedep=false, progressbar=true, onthefly=Dict(), kwargs...)
    A0=deepcopy(A)
    H0=deepcopy(H)
    data = Dict{String,Any}()

    numsteps = length(collect(0:dt:tmax))-1
    times = [(i-1)*dt for i=1:numsteps+1]
    
    @printf("prec : %.3e \n", prec)

    exp = measure(A0, obs; t=times[1])
    for i=1:length(obs)
        push!(data, obs[i].name => reshape(exp[i], size(exp[i])..., 1))
    end
    if reduceddensity
        :Nrho in keys(kwargs) ? Nrho = kwargs[:Nrho] : Nrho = 1
        exprho = rhoreduced_1site(A0,Nrho)
        push!(data, "Reduced ρ" => reshape(exprho, size(exprho)..., 1))
    end

    bonds = bonddims(A0)
    savebonddims && push!(data, "bonddims" => reshape([bonds...], length(bonds), 1))

    error && (errs = Vector{Float64}(undef, numsteps))
    timed && (ttdvp = Vector{Float64}(undef, numsteps))
    timed && (tproj = Vector{Float64}(undef, numsteps))
    effects && (efft = Vector{Any}(undef, numsteps))
    onthefly = copy(onthefly)
    ontheflyplotbool = !isempty(onthefly) && !isnothing(onthefly[:plot_obs])
    ontheflysavebool =  !isempty(onthefly) && !isempty(onthefly[:save_obs])
    ontheflysavebool && (onthefly[:save_obs] = intersect(onthefly[:save_obs], [ob.name for ob in obs]))

    F=nothing
    Afull=nothing
    iter = progressbar ? ProgressBar(numsteps; ETA=false) : 1:numsteps
    for tstep in iter
        if timedep
           Ndrive = kwargs[:Ndrive]
           Htime = kwargs[:Htime]
           if length(Ndrive)==1
              size(H0[Ndrive][1,1,:,:])==size(Htime[tstep][:,:]) ? nothing : throw(error("The size of Htime does not match the size of the non-interacting part of H at Ndrive"))
              H0[Ndrive][1,1,:,:] = H[Ndrive][1,1,:,:] + Htime[tstep][:,:]
           else
              for i=1:length(Ndrive)
                 site=Ndrive[i]
                 size(H0[site][end,1,:,:])==size(Htime[site][tstep][:,:]) ? nothing : throw(error("The size of Htime does not match the size of the non-interacting part of H at Ndrive=$site"))
                 H0[site][end,1,:,:] = H[site][end,1,:,:] + Htime[site][tstep][:,:]
              end
            end
        end

        A0, Afull, F, info = tdvp1sweep_dynamic!(dt, A0, H0, Afull, F;
                                                 obs=obs,
                                                 prec=prec,
                                                 Dlim=Dlim,
                                                 Dplusmax=Dplusmax,
                                                 timed=timed,
                                                 error=error,
                                                 kwargs...
                                                 )
        
        exp = info["obs"]
        bonds = info["dims"]
        maxbond = max(bonds...)
        progressbar && (iter.Dmax = maxbond)

        effects && (efft[tstep] = reduce(hcat, info["effect"]))
        error && (errs[tstep] = info["err"])
        timed && (ttdvp[tstep] = info["t2"] + info["t3"])
        timed && (tproj[tstep] = info["t1"])

        if tstep != 1
            for (i, ob) in enumerate(obs)
                data[ob.name] = cat(data[ob.name], exp[i]; dims=ndims(exp[i])+1)
            end
            if reduceddensity
                exprho = rhoreduced_1site(A0,Nrho)
            	data["Reduced ρ"] = cat(data["Reduced ρ"], exprho; dims=ndims(exprho)+1)
            end
            ontheflyplotbool && (tstep-1)%onthefly[:step] == 0 && ontheflyplot(onthefly, tstep-1, times, data)
            ontheflysavebool && (tstep-1)%onthefly[:step] == 0 && ontheflysave(onthefly, tstep-1, times, data)
        end
        if savebonddims
            data["bonddims"] = cat(data["bonddims"], bonds, dims=2)
        end
    end
    exp = measure(A0, obs; t=times[end])
    for (i, ob) in enumerate(obs)
        data[ob.name] = cat(data[ob.name], exp[i]; dims=ndims(exp[i])+1)
    end
    if reduceddensity
        exprho = rhoreduced_1site(A0,Nrho)
        data["Reduced ρ"] = cat(data["Reduced ρ"], exprho; dims=ndims(exprho)+1)
    end

    
    if effects
        efftarray = efft[1]
        for tstep=2:numsteps
            efftarray = cat(efftarray, efft[tstep], dims=3)
        end
        push!(data, "effects" => efftarray)
    end
    timed && push!(data, "projtime" => tproj)
    timed && push!(data, "tdvptime" => ttdvp)
    error && push!(data, "projerr" => errs)
    push!(data, "times" => times)
    return A0, data
end
