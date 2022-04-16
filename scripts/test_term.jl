using DrWatson
@quickactivate "Tesi"

includet(srcdir("CabibboMarinari.jl"))

using QCD, Plots, Statistics, CSV
using Term.progress

#=
when calling this function, the blue header shows up
but the progress bar doesn't. The bar is shown at 100%
completion after the loop finishes 
=#
function termalization!(L::Lattice)
	@track for n in 1:100
		lattice_overrelaxation!(L, 3)
		lattice_heatbath!(L, 1.0)
		if n % 10 == 0
			lattice_normalization!(L)
		end
	end
end

##

L = Lattice(4, 4, 4)
termalization!(L)