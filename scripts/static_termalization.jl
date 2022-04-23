using DrWatson
@quickactivate "Tesi"

using StaticQCD
nworkers() == 1 && initprocs(4)
@everywhere using StaticQCD

using Plots, Statistics, ProgressMeter, DataFrames, Measurements
import DelimitedFiles, CSV

include(scriptsdir("parameters.jl"))


#* ===== UTILITIES =====

function DrWatson._wsave(filename, data::Dict)
	if splitext(filename)[2] == ".dat" 
		DelimitedFiles.writedlm(filename, data, " = ")
	else
		save(filename, data)
	end
end

function DrWatson._wsave(filename, data::DataFrame)
	CSV.write(filename, data)
end

showall(x) = begin show(stdout, "text/plain", x); println() end

takeobservables(x::Pair) = (x.first,), (x.second,)
takeobservables(x) = first.(x), last.(x)

function incrementalmean(v, offset::Integer = 1)
	if offset ∉ 1:length(v)
		throw(ArgumentError("Given offset = $offset must be positive and less than or equal to length(v) = $(length(v))"))
	end

	[mean(v[offset:i]) for i in offset:length(v)], [std(v[offset:i]) for i in offset:length(v)]
end


#* ===== TERMALIZATION =====

function termalization!(L, params)
	@unpack β, nterm, nover, nnorm = params
	
	@showprogress 1 "Distributed Termalization" for n in 1:nterm
		one_termalization!(L, nover, β, n % nnorm == 0)
	end
end

function termalization(params)
	@unpack dims, latticestart = params
	L = Lattice(dims, start = latticestart)
	termalization!(L, params)
	L
end

function termalization!(L, params, observable::Function, v::Vector)
	@unpack β, nterm, nover, nnorm, nobs = params
	
	@showprogress 1 "Distributed Termalization" for n in 1:nterm
		one_termalization!(L, nover, β, n % nnorm == 0)
		n % nobs == 0 && push!(v, observable(L))
	end
	nterm % nobs ≠ 0 && push!(v, observable(L))
end

function termalization(params, observable::Function)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart)
	v = Float64[]
	termalization!(L, params, observable, v)
	v, L
end

function termalization!(L, params, observables, v)
	@unpack β, nterm, nover, nnorm, nobs = params
	
	@showprogress 1 "Distributed Termalization" for n in 1:nterm
		one_termalization!(L, nover, β, n % nnorm == 0)
		n % nobs == 0 && push!(v, [obs(L) for obs in observables])
	end

	# make sure the last iteration is measured
	nterm % nobs ≠ 0 && push!(v, [obs(L) for obs in observables])
end

function termalization(params, observables)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart)
	v = Vector{Float64}[]
	termalization!(L, params, observables, v)
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

		obsmeasurements, = termalization(params, obsfunctions)

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

	sort!(df, [:dims, :β])

	if save_df
		dfname = savename(allparams, "csv")
		safesave(datadir(folder, dfname), df)
	end

	df
end

run() = begin include(scriptsdir("parameters.jl")); run(TermParams(), ObsParams()) end
run(s::String) = begin include(scriptsdir("parameters.jl")); run(TermParams(), ObsParams(), s) end