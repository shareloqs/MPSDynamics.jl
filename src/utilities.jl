mutable struct ProgressBar
    numsteps::Int
    ETA::Bool
    times::Vector{<:Float64}
    Dmax::Int
    length::UInt
end

"""
    ProgressBar(numsteps::Int; ETA=false, last=10)

An iterable returning values from 1 to `numsteps`. Displays a progress bar for the for loop where it has been called.
If ETA is true then displays an estimation of the remaining time calculated based on the time spent computing the last `last` values.
The progress bar can be activated or deactivated by setting the `progressbar` keyword argument in `runsim` to true or false."""
ProgressBar(numsteps::Int; ETA=false, last=10) = ProgressBar(numsteps, ETA, fill(0., (ETA ? last : 1)+1), 0, displaysize(stdout)[2] > 54 ? min(displaysize(stdout)[2]-54, 50) : error("Error : Terminal window too narrow"))

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
        dtETA = (tnow - bar.times[2+state%Ntimes])*(bar.numsteps - state)/min(state-1, Ntimes)
        dtiter = tnow - bar.times[2+(state-1)%Ntimes]
        bar.times[2+state%Ntimes] = tnow
        elapsedstr = Dates.format(Time(0)+Second(floor(Int, dtelapsed)), dtelapsed>3600 ? "HH:MM:SS" : "MM:SS")
        ETAstr = Dates.format(Time(0)+Second(floor(Int, dtETA)), dtETA>3600 ? "HH:MM:SS" : "MM:SS")
        iterstr = Dates.format(Time(0)+Millisecond(floor(Int, 1000*dtiter)), dtiter>60 ? "MM:SS.sss" : "SS.sss")
        print("\r")
        printstyled("$(round(100*state/bar.numsteps, digits=1))% "; color = :green, bold=true)
        print("┣"*"#"^(round(Int, state/bar.numsteps*bar.length))*" "^(bar.length-round(Int, state/bar.numsteps*bar.length))*"┫")
        printstyled(" $state/$(bar.numsteps) [$(elapsedstr)s"*(bar.ETA ? "; ETA:$(ETAstr)s" : "")*"; $(iterstr)s/it"*(bar.Dmax > 0 ? "; Dmax=$(bar.Dmax)" : "")*"]"; color = :green, bold=true)
    end
    return (state, state+1)
end

"""
    onthefly(;plot_obs=nothing::Union{<:Observable, Nothing}, save_obs=Vector{Observable}(undef, 0)::Union{<:Observable, Vector{<:Observable}}, savedir="auto", step=10::Int, func=identity<:Function, compare=nothing::Union{Tuple{Vector{Float64}, Vector{Float64}}, Nothing}, clear=identity<:Function)

Helper function returning a dictionnary containing the necessary arguments for on the fly plotting or saving in the `runsim` function.

# Arguments

* `plot_obs` : Observable to plot
* `save_obs` : List of Observable(s) to save
* `savedir` : Used to specify the path where temporary files are stored, default is `"auto"` which saves in a "tmp" folder in the run folder (generally located at "~/MPSDynamics/<unid>/tmp/"). 
* `step` : Number of time steps every which the function plots or saves the data
* `func` : Function to apply to the result of measurement of plot_obs (for example `real` or `abs` to make a complex result possible to plot)
* `compare` : Tuple `(times, data)` of previous results to compare against the plot_obs results,
* `clear` : Function used to clear the output before each attempt to plot, helpful for working in Jupyter notebooks where clear=IJulia.clear_output allows to reduce the size of the cell output

# Examples

For example in the [Spin-Boson model example](@ref "The Spin Boson Model"), adding the following argument to the `runsim` function allows to save
the "sz" observable `ob1` in the directory "~/MPSDynamics/<unid>/tmp/" and plot its real part during the simulation:
```julia
runsim(..., onthefly=onthefly(plot_obs=ob1, save_obs=[ob1], savedir="auto", step=10, func=real))
```
To merge the temporary files in one usable file, one can then use [`MPSDynamics.mergetmp`](@ref).
"""
function onthefly(;plot_obs=nothing::Union{<:Observable, Nothing}, save_obs=Vector{Observable}(undef, 0)::Union{<:Observable, Vector{<:Observable}}, savedir="auto", step=10::Int, func=identity::Function, compare=nothing::Union{Tuple{Vector{Float64}, Vector{Float64}}, Nothing}, clear=nothing)
    if isnothing(plot_obs) && isempty(save_obs)
        error("Must provide an observable to plot/save")
    end
    # !isempty(save_obs) && (isdir(savedir) || mkdir(savedir))
    plt = isnothing(plot_obs) ? nothing : plot(title="Intermediate Results", xlabel="t", ylabel=plot_obs.name)
    !isnothing(compare) && plot!(compare)
    println("On the fly mode activated")
    return Dict(:plot_obs => isnothing(plot_obs) ? nothing : plot_obs.name, :save_obs => [ob.name for ob in save_obs], :savedir => savedir, :step => step, :func => func, :clear => clear, :compare => compare, :plot => plt)
end

"""
    ontheflyplot(onthefly, tstep, times, data)

Plots data according to the arguments of the `onthefly` dictionnary."""
function ontheflyplot(onthefly, tstep, times, data)
    tstep < onthefly[:step]+2 && plot!(onthefly[:plot], [], [])
    N = ndims(data[onthefly[:plot_obs]])
    N == 1 ? slicefunc(arr, N, i) = arr[i] : slicefunc = selectdim
    append!(onthefly[:plot].series_list[end], times[tstep-onthefly[:step]+1:tstep+1], map(onthefly[:func], (slicefunc(data[onthefly[:plot_obs]], N, i) for i in (tstep-onthefly[:step]+1):(tstep+1))))
    !isnothing(onthefly[:clear]) && onthefly[:clear](true)
    sleep(0.05);display(onthefly[:plot])
end

"""
    ontheflysave(onthefly, tstep, times, data)

Saves data according to the arguments of the `onthefly` dictionnary."""
function ontheflysave(onthefly, tstep, times, data)
    jldopen(onthefly[:savedir]*"tmp$(tstep÷onthefly[:step]).jld", "w") do file
        write(file, "times", times[tstep-onthefly[:step]+1:tstep+1])
        for name in onthefly[:save_obs]
            write(file, name, data[name][tstep-onthefly[:step]+1:tstep+1])
        end
    end
end

"""
    mergetmp(tmpdir; fields=[], overwrite=true)

Merges the temporary files created by the `ontheflysave` function at the directory `tmpdir` and returns a dictionnary containing the resulting data.
By default all fields are present but one can select the fields of interest with a list of names in `fields`."""
function mergetmp(tmpdir; fields=[], overwrite=true)
    tmpdir[end] != '/' && (tmpdir *= '/')
    !isdir(tmpdir) && error("Choose a valid directory")
    files = [walkdir(tmpdir)...][1][3]
    files = filter(x->!isnothing(match(r"tmp(\d+).jld", x)), files)
    files = sort(files, by=(x -> parse(Int, match(r"(\d+)", x).captures[1])))
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