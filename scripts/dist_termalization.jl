using DrWatson
@quickactivate "Tesi"

using Distributed
addprocs(4, exeflags="--project")

using QCD, Plots, Statistics, ProgressMeter, DataFrames
import DelimitedFiles, CSV

##

@everywhere begin
	using DrWatson
	include(srcdir("CabibboMarinari.jl"))
	include(scriptsdir("parameters.jl"))
end

##

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


#* ===== TERMALIZATION =====

function one_termalization!(L::Lattice, nover::Integer, β::Real, normalize::Bool = false)
	dist_lattice_overrelaxation!(L, nover)
	dist_lattice_heatbath!(L, β)
	normalize && lattice_normalization!(L)
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
	@unpack β, nterm, nover, nnorm = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		one_termalization!(L, nover, β, n % nnorm == 0)
		push!(v, observable(L))
	end
end

function termalization(params, observable::Function)
	@unpack dims, latticestart, sp2type = params
	lattice = Lattice(dims..., start = latticestart, type = sp2type)
	v = []
	termalization!(lattice, params, observable, v)
	v, lattice
end

function termalization!(L::Lattice, params, observables, v)
	@unpack β, nterm, nover, nnorm = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		one_termalization!(L, nover, β, n % nnorm == 0)

		w = []
		for obs in observables
			push!(w, obs(L))
		end
		push!(v, w)
	end
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
		@unpack meanoffset = params

		display(params)
		d = Dict(params)

		measurements, = termalization(params, obsfunctions)

		for (measurement, obsname) in zip(eachcol(measurements), obsnames)
			obsmean, xrange = incrementalmean(measurement, meanoffset)
			d[obsname] = obsmean[end] # add final mean to dictionary

			println("mean $obsname = $(obsmean[end])")

			if display_plot || save_plot
				plottitle = savename(params, connector = ", ", sort = false)
				p = plot(measurement, label = obsname, title = plottitle, titlefontsize = 10)
				plot!(p, xrange, obsmean, label = "mean $obsname")
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