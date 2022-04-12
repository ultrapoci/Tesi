using DrWatson
@quickactivate "Tesi"

includet(srcdir("CabibboMarinari.jl"))

using QCD, ProgressMeter, Plots

##

function termalization!(L::Lattice, avg_plaqs::Vector{Float64}, params::Dict)
	@unpack β, nterm, nover, norm_every = params
	
	@showprogress 1 "Termalization" for n in 1:nterm
		#lattice_overrelaxation!(L, nover)
		lattice_heatbath!(L, β)
		if n % norm_every == 0
			lattice_normalization!(L)
		end
		push!(avg_plaqs, averageplaquette(L))
	end
end


function termalization(params::Dict)
	@unpack dims, lattice_start = params
	lattice = Lattice(dims, start = lattice_start)
	avg_plaqs = Float64[]
	termalization!(lattice, avg_plaqs, params)
	lattice, avg_plaqs
end

##

include("parameters.jl")
display(params)
lattice, avg_plaqs = termalization(params);
plot(avg_plaqs)
