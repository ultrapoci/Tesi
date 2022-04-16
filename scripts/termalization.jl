using DrWatson
@quickactivate "Tesi"

includet(srcdir("CabibboMarinari.jl"))

using QCD, Plots, Statistics, CSV, ProgressMeter

function termalization!(L::Lattice, params::Dict)
	@unpack β, nterm, nover, nnorm = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		lattice_overrelaxation!(L, nover)
		lattice_heatbath!(L, β)
		if n % nnorm == 0
			lattice_normalization!(L)
		end
	end
end

function termalization!(L::Lattice, params::Dict, observable::Function, v::Vector)
	@unpack β, nterm, nover, nnorm = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		lattice_overrelaxation!(L, nover)
		lattice_heatbath!(L, β)
		if n % nnorm == 0
			lattice_normalization!(L)
		end
		push!(v, observable(L))
	end
end

function termalization!(L::Lattice, params::Dict, observables, v)
	@unpack β, nterm, nover, nnorm = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
	for n in 1:nterm
		lattice_overrelaxation!(L, nover)
		lattice_heatbath!(L, β)
		if n % nnorm == 0
			lattice_normalization!(L)
		end

		w = []
		for obs in observables
			push!(w, obs(L))
		end
		push!(v, w)
	end
end

function termalization(params::Dict)
	@unpack dims, latticestart = params
	lattice = Lattice(dims, start = lattice_start)
	termalization!(lattice, params)
	lattice
end

function termalization(params::Dict, observables)
	@unpack dims, latticestart, sp2type = params
	lattice = Lattice(dims..., type = sp2type, start = latticestart)
	v = []
	termalization!(lattice, params, observables, v)
	lattice, reduce(hcat, v)
end

function incrementalmean(v)
	res = []
	for i in 1:length(v)
		push!(res, mean(v[1:i]))
	end
	res
end

##

include("parameters.jl")

for params in dict_list(allparams)
	show(stdout, "text/plain", params)
	println()

	@unpack observables, to_plot = params

	lattice, all_obs = termalization(params, observables)
	d = copy(params)
	for (obs, observable) in zip(eachrow(all_obs), observables)
		obs_mean = incrementalmean(obs)
		d[Symbol(observable)] = obs_mean[end]
		obs_name = String(Symbol(observable))

		println("mean $obs_name = $(obs_mean[end])")

		if to_plot
			plottitle = savename(params, ignores = "to_plot")
			p = plot(obs, label = obs_name, title = plottitle, titlefontsize = 10)
			plot!(p, obs_mean, label = "mean $obs_name")
			plotname = savename(obs_name, params, "png", ignores = "to_plot")		
			savefig(p, plotsdir(plotname))
			display(p)
		end
	end

	delete!(d, :to_plot)
	delete!(d, :observables)
	csvname = savename(params, "csv", ignores = "to_plot")
	CSV.write(datadir(csvname), d, writeheader = false)
	println()
end