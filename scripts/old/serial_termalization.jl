using DrWatson
@quickactivate "Tesi"

using QCD, Plots, Statistics, ProgressMeter, DataFrames, Measurements

include(srcdir("Utilities.jl"))
include(scriptsdir("parameters.jl"))


#* ===== TERMALIZATION =====

function one_termalization!(L::Lattice, nover::Integer, β::Real, normalize::Bool = false; log = false)
	lattice_overrelaxation!(L, nover, log = log)
	lattice_heatbath!(L, β, log = log)
	normalize && lattice_normalization!(L, log = log)
end

function termalization!(L::Lattice, params)
	@unpack β, nterm, nover, nnorm = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		one_termalization!(L, nover, β, n % nnorm == 0)
	end
end

function termalization(params)
	@unpack dims, latticestart, sp2type = params
	lattice = Lattice(dims, start = latticestart, type = sp2type)
	termalization!(lattice, params)
	lattice
end

function termalization!(L::Lattice, params, observable::Function, v::Vector)
	@unpack β, nterm, nover, nnorm, nobs = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		one_termalization!(L, nover, β, n % nnorm == 0)
		n % nobs == 0 && push!(v, observable(L))
	end
	nterm % nobs ≠ 0 && push!(v, observable(L))
end

function termalization(params, observable::Function)
	@unpack dims, latticestart, sp2type = params
	lattice = Lattice(dims..., start = latticestart, type = sp2type)
	v = []
	termalization!(lattice, params, observable, v)
	v, lattice
end

function termalization!(L::Lattice, params, observables, v)
	@unpack β, nterm, nover, nnorm, nobs = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		one_termalization!(L, nover, β, n % nnorm == 0)

		n % nobs == 0 && push!(v, [obs(L) for obs in observables])
	end

	nterm % nobs ≠ 0 && push!(v, [obs(L) for obs in observables])
end

function termalization(params, observables)
	@unpack dims, latticestart, sp2type = params
	lattice = Lattice(dims..., type = sp2type, start = latticestart)
	v = []
	termalization!(lattice, params, observables, v)
	reduce(hcat, v)', lattice
end


#* ===== RUN =====

function run(allparams::TermParams, obsparams::ObsParams, folder = "")
	@unpack observables, save_plot, display_plot, save_dat, save_jld2, save_df = obsparams
	obsnames, obsfunctions = takeobservables(observables)

	df = DataFrame()

	for params in dict_list(allparams)
		@unpack startobs = params

		display(params)
		d = Dict(params)

		obsmeasurements, = termalization(params, obsfunctions)

		for (obsmeasurement, obsname) in zip(eachcol(obsmeasurements), obsnames)
			obsmean, obserror = incremental_measurement(obsmeasurement, startobs)
			d[obsname] = measurement(obsmean[end], obserror[end]) # add final mean and std to dictionary

			println("$obsname = $(d[obsname])")

			if display_plot || save_plot
				plottitle = savename(params, connector = ", ", sort = false)
				p = plot(obsmeasurement, label = obsname, title = plottitle, titlefontsize = 10);
				plot!(p, startobs:length(obsmeasurement), obsmean, label = "mean $obsname");
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

	sort!(df, [:dims, :β])

	if save_df
		dfname = savename(allparams, "csv")
		safesave(datadir(folder, dfname), df)
	end

	df
end

run() = begin include(scriptsdir("parameters.jl")); run(TermParams(), ObsParams()) end
run(s::String) = begin include(scriptsdir("parameters.jl")); run(TermParams(), ObsParams(), s) end