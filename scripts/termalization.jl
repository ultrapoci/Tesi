using DrWatson
@quickactivate "Tesi"

if "TERMALIZATION_LOG" ∉ keys(ENV)
	ENV["TERMALIZATION_LOG"] = "0"
end

if "TERMALIZATION_NPROCS" ∉ keys(ENV)
	ENV["TERMALIZATION_NPROCS"] = "4"
end

using DistributedQCD
nworkers() == 1 && initprocs(parse(Int, ENV["TERMALIZATION_NPROCS"]))
@everywhere using DrWatson 
@everywhere @quickactivate "Tesi"
@everywhere using DistributedQCD

using Plots, DataFrames, Measurements

include(srcdir("Utilities.jl"))
include(scriptsdir("parameters.jl"))


#* ===== TERMALIZATION =====

function termalization!(L, params; log = false)
	@unpack β, nterm, nover, nnorm = params

	pbar = getpbar(nterm, desc = " Distributed Termalization", enabled = !log)

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		next!(pbar, showvalues = generate_showvalues(:iter => n, :total => nterm))
	end
end

function termalization(params; kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart)
	termalization!(L, params; kwargs...)
	L
end

function termalization!(L, params, observable::Function, v::Vector; log = false)
	@unpack β, nterm, nover, nnorm, nobs = params

	pbar = getpbar(nterm, desc = " Distributed Termalization", enabled = !log)

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		n % nobs == 0 && push!(v, observable(L; log = log, iter = n))
		next!(pbar, showvalues = generate_showvalues(:iter => n, :total => nterm))
	end
	
	# make sure the last iteration is measured
	nterm % nobs ≠ 0 && push!(v, observable(L; log = log, iter = n))
end

function termalization(params, observable::Function; kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart)
	v = Float64[]
	termalization!(L, params, observable, v; kwargs...)
	v, L
end

function termalization!(L, params, observables, v; log = false)
	@unpack β, nterm, nover, nnorm, nobs = params

	pbar = getpbar(nterm, desc = " Distributed Termalization", enabled = !log)

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		n % nobs == 0 && push!(v, [obs(L; log = log, iter = n) for obs in observables])
		next!(pbar, showvalues = generate_showvalues(:iter => n, :total => nterm))
	end

	# make sure the last iteration is measured
	nterm % nobs ≠ 0 && push!(v, [obs(L; log = log, iter = n) for obs in observables])
end

function termalization(params, observables; kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart)
	v = Vector{Float64}[]
	termalization!(L, params, observables, v; kwargs...)
	reduce(hcat, v)', L
end


#* ===== RUN =====

function run(allparams::TermParams, obsparams::ObsParams, folder = "")
	@unpack observables, save_plot, display_plot, save_dat, save_jld2, save_df = obsparams
	obsnames, obsfunctions = takeobservables(observables)

	df = DataFrame()

	for params in dict_list(allparams)
		@unpack meanoffset = params

		display(params)
		d = Dict(params)

		obsmeasurements, = termalization(params, obsfunctions, log = ENV["TERMALIZATION_LOG"] == "1")

		for (obsmeasurement, obsname) in zip(eachcol(obsmeasurements), obsnames)
			obsmean, obserror = incrementalmean(obsmeasurement, meanoffset)
			d[obsname] = measurement(obsmean[end], obserror[end]) # add final mean and std to dictionary

			println("$obsname = $(d[obsname])")

			if display_plot || save_plot
				plottitle = savename(params, connector = ", ", sort = false)
				p = plot(obsmeasurement, label = obsname, title = plottitle, titlefontsize = 10);
				plot!(p, meanoffset:length(obsmeasurement), obsmean, label = "mean $obsname");
				plotname = savename(obsname, params, "png", sort = false)	
				save_plot && safesave(plotsdir(folder, plotname), p)
				display_plot && display(p)
			end
		end

		if save_jld2
			jld2name = savename(params, "jld2", sort = false)
			safesave(datadir(folder, "jld2", jld2name), d)
		end 

		if save_dat
			datname = savename(params, "dat", sort = false)
			safesave(datadir(folder, "dat", datname), d)
		end

		append!(df, d)
		println()
	end

	# remove columns of missing values
	df = df[!, (x->eltype(x)!=Missing).(eachcol(df))]
	sort!(df, [:dims, :β])

	if save_df
		dfname = savename(allparams, "csv")
		safesave(datadir(folder, dfname), df)
	end

	df
end

run() = begin include(scriptsdir("parameters.jl")); run(TermParams(), ObsParams()) end
run(s::String) = begin include(scriptsdir("parameters.jl")); run(TermParams(), ObsParams(), s) end