mutable struct ProgressBar
    numsteps::Int
    ETA::Bool
    times::Vector{<:Float64}
    Dmax::Int
end

"""An iterable returning values from 1 to `numsteps`. Displays a progress bar for the for loop where it has been called. If ETA is true then displays an estimation of the remaining time calculated based on the time spent computing the last `last` values."""
ProgressBar(numsteps::Int; ETA=false, last=10) = ProgressBar(numsteps, ETA, fill(0., (ETA ? last : 1)+1), 0)

function Base.iterate(bar::ProgressBar, state=1)
    Ntimes = length(bar.times)-1
    if state > bar.numsteps
        println()
        return nothing
    elseif state == 1
        printstyled("\nCompiling..."; color=:red, bold=true)
        bar.times = fill(time(), Ntimes+1)
    else
        tnow = time()
        dtelapsed = tnow - bar.times[1]
        dtETA = tnow - bar.times[2+state%Ntimes]
        dtiter = tnow - bar.times[2+(state-1)%Ntimes]
        bar.times[2+state%Ntimes] = tnow
        elapsedstr = Dates.format(Time(0)+Second(floor(Int, dtelapsed)), dtelapsed>3600 ? "HH:MM:SS" : "MM:SS")
        ETAstr = Dates.format(Time(0)+Second(floor(Int, dtETA*(bar.numsteps - state)/min(state-1, Ntimes))), dtETA>3600 ? "HH:MM:SS" : "MM:SS")
        iterstr = Dates.format(Time(0)+Millisecond(floor(Int, 1000*dtiter)), dtiter>60 ? "MM:SS.sss" : "SS.sss")

        printstyled("\r$(round(100*state/bar.numsteps, digits=1))% "; color = :green, bold=true)
        print("┣"*"#"^(round(Int, state/bar.numsteps*50))*" "^(50-round(Int, state/bar.numsteps*50))*"┫")
        printstyled(" $state/$(bar.numsteps) [$(elapsedstr)s"*(bar.ETA ? "; ETA:$(ETAstr)s" : "")*"; $(iterstr)s/it"*(bar.Dmax > 0 ? "; Dmax=$(bar.Dmax)" : "")*"]"; color = :green, bold=true)
    end
    return (state, state+1)
end

"""
    onthefly(;plot_obs=nothing::Union{<:Observable, Nothing}, save_obs=Vector{Observable}(undef, 0)::Union{<:Observable, Vector{<:Observable}}, savedir=string(homedir(),"/MPSDynamics/tmp/"), step=10::Int, func=identity<:Function, compare=nothing::Union{Tuple{Vector{Float64}, Vector{Float64}}, Nothing}, clear=identity<:Function)

Helper function to send arguments for on the fly plotting or saving.
    `plot_obs` : Observable to plot
    `save_obs` : Observable(s) to save
    `savedir` : Used to specify the path where temporary files are stored
    `step` : Number of time steps every which the function plots or saves the data
    `func` : Function to apply to the result of measurement of plot_obs (for example `real` or `abs` to make a complex result possible to plot)
    `compare` : Tuple `(times, data)` of previous results to compare against the plot_obs results,
    `clear` : Function used to clear the output before each attempt to plot, helpful for working in Jupyter notebooks where clear=IJulia.clear_output allows to reduce the size of the cell output
"""
function onthefly(;plot_obs=nothing::Union{<:Observable, Nothing}, save_obs=Vector{Observable}(undef, 0)::Union{<:Observable, Vector{<:Observable}}, savedir=string(homedir(),"/MPSDynamics/tmp/"), step=10::Int, func=identity<:Function, compare=nothing::Union{Tuple{Vector{Float64}, Vector{Float64}}, Nothing}, clear=nothing)
    if isnothing(plot_obs) && isempty(save_obs)
        error("Must provide an observable to plot/save")
    end
    !isempty(save_obs) && (isdir(savedir) || mkdir(savedir))
    plt = isnothing(plot_obs) ? nothing : plot(title="Intermediate Results", xlabel="t", ylabel=plot_obs.name)
    println("On the fly mode activated")
    return Dict(:plot_obs => plot_obs.name, :save_obs => [ob.name for ob in save_obs], :savedir => savedir, :step => step, :func => func, :clear => clear, :compare => compare, :plot => plt)
end

"""
    ontheflyplot(onthefly, tstep, times, data)

Plots data according to the arguments of the `onthefly` dictionnary"""
function ontheflyplot(onthefly, tstep, times, data)
    tstep < onthefly[:step]+2 && plot!(onthefly[:plot], [], [])
    append!(onthefly[:plot].series_list[end], times[tstep-onthefly[:step]+1:tstep+1], onthefly[:func].(data[onthefly[:plot_obs]][tstep-onthefly[:step]+1:tstep+1]))
    !isnothing(onthefly[:compare]) && plot!(onthefly[:compare])
    !isnothing(onthefly[:clear]) && onthefly[:clear](true)
    sleep(0.05);display(onthefly[:plot])
end

"""
    ontheflysave(onthefly, tstep, times, data)

Saves data according to the arguments of the `onthefly` dictionnary"""
function ontheflysave(onthefly, tstep, times, data)
    jldopen(onthefly[:savedir]*"tmp$(tstep÷onthefly[:step]).jld", "w") do file
        write(file, "times", times[tstep-onthefly[:step]+1:tstep+1])
        for name in onthefly[:save_obs]
            write(file, name, data[name][tstep-onthefly[:step]+1:tstep+1])
        end
    end
end

"""
    mergetmp(;tmpdir=string(homedir(),"/MPSDynamics/tmp/"), fields=[], overwrite=true)

Merges the temporary files created by the `ontheflysave` function and returns a dictionnary containing the resulting data"""
function mergetmp(;tmpdir=string(homedir(),"/MPSDynamics/tmp/"), fields=[], overwrite=true)
    tmpdir[end] != '/' && (tmpdir *= '/')
    !isdir(tmpdir) && error("Choose a valid directory")

    files = sort([walkdir(tmpdir)...][1][3], by=(x->parse(Int, match(r"(\d+)", x).captures[1])))
    isempty(fields) && (fields = keys(JLD.load(tmpdir*files[1])))
    merged_data = Dict(ob => [] for ob in fields)
    i = 1
    for file in files
        for ob in fields
            append!(merged_data[ob], JLD.load(tmpdir*file, ob))
        end
        i += 1
    end
    if overwrite
        for file in files
            rm(tmpdir*file)
        end
        JLD.save(tmpdir*files[end], Iterators.flatten(merged_data)...)
    end
    return merged_data
end