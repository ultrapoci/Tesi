using DrWatson
@quickactivate "Tesi"

using QCD, Plots, Statistics, ProgressMeter
import DelimitedFiles
includet(srcdir("CabibboMarinari.jl"))

function DrWatson._wsave(filename, data::Dict)
	if splitext(filename)[2] == ".dat" 
		DelimitedFiles.writedlm(filename, data, " = ")
	else
		save(filename, data)
	end
end

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

function incrementalmean(v, offset::Integer = 1)
	if offset ∉ 1:length(v)
		throw(ArgumentError("Given offset = $offset must be positive and less than or equal to length(v) = $(length(v))"))
	end

	[mean(v[offset:i]) for i in offset:length(v)], offset:length(v)
end

showall(x) = begin show(stdout, "text/plain", x); println() end

iterobs(ob::Function) = (String(Symbol(ob)),)
iterobs(ob::Tuple) = String.(Symbol.(ob))

##

include("parameters.jl")
@unpack observables, to_plot, meanoffset = obsparams

for params in dict_list(allparams)
	display(params)
	d = copy(params)
	d["observables"] = Dict()

	measurements, = termalization(params, observables)

	for (measurement, obs_name) in zip(eachcol(measurements), iterobs(observables))
		obs_mean, xrange = incrementalmean(measurement, meanoffset)
		d["observables"][obs_name] = obs_mean[end] # add final mean to dictionary

		println("mean $obs_name = $(obs_mean[end])")

		if to_plot
			plottitle = savename(params)
			p = plot(measurement, label = obs_name, title = plottitle, titlefontsize = 10)
			plot!(p, xrange, obs_mean, label = "mean $obs_name")
			plotname = savename(obs_name, params, "png")	
			savefig(p, plotsdir(plotname))
			display(p)
		end
	end

	jld2name = savename(params, "jld2")
	wsave(datadir("jld2", jld2name), d)

	datname = savename(params, "dat")
	wsave(datadir("dat", datname), merge(params, d["observables"]))

	println()
end