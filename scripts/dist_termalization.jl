using DrWatson
@quickactivate "Tesi"

using Distributed
nworkers() == 1 && addprocs(3, exeflags="--project")

##

@everywhere using DrWatson

@everywhere begin
	using QCD, Plots, Statistics, ProgressMeter, DataFrames
	import DelimitedFiles, CSV

	include(srcdir("CabibboMarinari.jl"))
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

		[mean(v[offset:i]) for i in offset:length(v)], offset:length(v)
	end

	function nsteps(params)
		@unpack nterm = params
		n = dict_list_count(params) ÷ length(nterm)
		sum(n .* nterm)
	end

	#* ===== TERMALIZATION =====

	function one_termalization!(L::Lattice, nover::Integer, β::Real, normalize::Bool = false)
		lattice_overrelaxation!(L, nover)
		lattice_heatbath!(L, β)
		normalize && lattice_normalization!(L)
	end

	function termalization!(L::Lattice, params)
		@unpack β, nterm, nover, nnorm = params
		
		for n in 1:nterm
			one_termalization!(L, nover, β, n % nnorm == 0)
		end
	end

	function termalization(params)
		@unpack dims, latticestart, sp2type = params
		lattice = Lattice(dims, start = latticestart, type = sp2type)
		termalization!(lattice, params)
		lattice
	end

	function termalization!(L::Lattice, params, observables, v)
		@unpack β, nterm, nover, nnorm, nobs = params
		
		for n in 1:nterm
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
end

#* ===== RUN =====

function dist_run(allparams::TermParams, obsparams::ObsParams, folder = "")
	@unpack observables, save_plot, display_plot, save_dat, save_jld2, save_df = obsparams
	obsnames, obsfunctions = takeobservables(observables)
	niter = try
		dict_list_count(allparams)
	catch
		1
	end
	v = pmap(enumerate(dict_list(allparams))) do (i, params)
		@info "Iteration $i/$niter"

		@unpack meanoffset = params

		d = Dict(params)

		measurements, = termalization(params, obsfunctions)

		for (measurement, obsname) in zip(eachcol(measurements), obsnames)
			obsmean, xrange = incrementalmean(measurement, meanoffset)

			d[obsname] = obsmean[end] # add final mean to dictionary

			println("mean $obsname = $(obsmean[end])")

			if display_plot || save_plot
				plottitle = savename(params, connector = ", ", sort = false)
				p = plot(measurement, label = obsname, title = plottitle, titlefontsize = 10);
				plot!(p, xrange, obsmean, label = "mean $obsname");
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

		d
	end

	df = DataFrame()
	for d in v
		append!(df, d)
	end

	sort!(df, [:dims, :β])

	if save_df
		dfname = savename(allparams, "csv")
		safesave(datadir(folder, dfname), df)
	end

	df
end

dist_run() = begin include(scriptsdir("parameters.jl")); run(TermParams(), ObsParams()) end