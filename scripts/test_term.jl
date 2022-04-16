using DrWatson
@quickactivate "Tesi"

includet(srcdir("CabibboMarinari.jl"))

# QCD is a local package
using QCD
using Term.progress
import ProgressMeter

#=
when calling this function, the blue header shows up
but the progress bar doesn't. The bar is shown at 100%
completion after the loop finishes.
=#
function test1!(L::Lattice)
	@track for n in 1:100
		lattice_overrelaxation!(L, 3)
		lattice_heatbath!(L, 1.0)
		if n % 10 == 0
			lattice_normalization!(L)
		end
	end
end

#=
ProgressMeter.jl works as intended, and shows that the loop
doesn't finish in so little time that the progress bar cannot
appear fast enough. 
=#
function test2!(L::Lattice)
	ProgressMeter.@showprogress 1 for n in 1:100
		lattice_overrelaxation!(L, 3)
		lattice_heatbath!(L, 1.0)
		if n % 10 == 0
			lattice_normalization!(L)
		end
	end
end

##

L = Lattice(4, 4, 4)
test1!(L)
test2!(L)