using DrWatson
@quickactivate "Tesi"

includet(srcdir("CabibboMarinari.jl"))

using QCD, ProgressMeter, Plots

##

function termalization!(L::Lattice, params::Dict)
	@unpack β, nterm, nover, norm_every = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		lattice_overrelaxation!(L, nover)
		# lattice_heatbath!(L, β)
		if n % norm_every == 0
			lattice_normalization!(L)
		end
	end
end

function termalization!(L::Lattice, params::Dict, observable::Function, v::Vector{Real})
	@unpack β, nterm, nover, norm_every = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		lattice_overrelaxation!(L, nover)
		# lattice_heatbath!(L, β)
		if n % norm_every == 0
			lattice_normalization!(L)
		end
		push!(v, observable(L))
	end
end


function termalization(params::Dict)
	@unpack dims, lattice_start = params
	lattice = Lattice(dims, start = lattice_start)
	termalization!(lattice, params)
	lattice
end

function termalization(params::Dict, observable::Function)
	@unpack dims, lattice_start, sp2type = params
	lattice = Lattice(dims..., type = sp2type, start = lattice_start)
	v = Real[]
	termalization!(lattice, params, observable, v)
	lattice, v
end

##

include("parameters.jl")
display(params)
lattice, v = termalization(params, averageplaquette);
println(averageplaquette(lattice))
plot(v)
