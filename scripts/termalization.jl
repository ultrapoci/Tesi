using DrWatson
@quickactivate "Tesi"

using QCD, Plots, Statistics, ProgressMeter
import DelimitedFiles
include(srcdir("CabibboMarinari.jl"))


#* ===== UTILITIES =====

function DrWatson._wsave(filename, data::Dict)
	if splitext(filename)[2] == ".dat" 
		DelimitedFiles.writedlm(filename, data, " = ")
	else
		save(filename, data)
	end
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
	lattice_overrelaxation!(L, nover)
	lattice_heatbath!(L, β/2) # β/2 is to have 8/g² as coefficient
	normalize && lattice_normalization!(L)
end

function termalization!(L::Lattice, params::Dict)
	@unpack β, nterm, nover, nnorm = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		one_termalization!(L, nover, β, n % nnorm == 0)
	end
end

function termalization(params::Dict)
	@unpack dims, latticestart, sp2type = params
	lattice = Lattice(dims, start = latticestart, type = sp2type)
	termalization!(lattice, params)
	lattice
end

function termalization!(L::Lattice, params::Dict, observable::Function, v::Vector)
	@unpack β, nterm, nover, nnorm = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		one_termalization!(L, nover, β, n % nnorm == 0)
		push!(v, observable(L))
	end
end

function termalization(params::Dict, observable::Function)
	@unpack dims, latticestart, sp2type = params
	lattice = Lattice(dims..., start = latticestart, type = sp2type)
	v = []
	termalization!(lattice, params, observable, v)
	v, lattice
end

function termalization!(L::Lattice, params::Dict, observables, v)
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

function termalization(params::Dict, observables)
	@unpack dims, latticestart, sp2type = params
	lattice = Lattice(dims..., type = sp2type, start = latticestart)
	v = []
	termalization!(lattice, params, observables, v)
	reduce(hcat, v)', lattice
end


#* ===== RUN =====

function run()
	include(scriptsdir("parameters.jl"))

	@unpack observables, to_plot, meanoffset = obsparams
	obsnames, obsfunctions = takeobservables(observables)

	for params in dict_list(allparams)
		display(params)
		d = copy(params)

		measurements, = termalization(params, obsfunctions)

		for (measurement, obsname) in zip(eachcol(measurements), obsnames)
			obsmean, xrange = incrementalmean(measurement, meanoffset)
			d[obsname] = obsmean[end] # add final mean to dictionary

			println("mean $obsname = $(obsmean[end])")

			if to_plot
				plottitle = savename(params, connector = ", ")
				p = plot(measurement, label = obsname, title = plottitle, titlefontsize = 10)
				plot!(p, xrange, obsmean, label = "mean $obsname")
				plotname = savename(obsname, params, "png", ignores = "latticestart")	
				safesave(plotsdir(plotname), p)
				display(p)
			end
		end

		jld2name = savename(params, "jld2", ignores = "latticestart")
		safesave(datadir("jld2", jld2name), d)

		datname = savename(params, "dat", ignores = "latticestart")
		safesave(datadir("dat", datname), d)

		println()
	end
end