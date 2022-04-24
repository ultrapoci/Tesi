Base.@kwdef struct ObsParams <: Params
	observables = "avg_plaq" => averageplaquette
	save_plot = true
	display_plot = true
	save_dat = true
	save_jld2 = true
	save_df = true
end

Base.@kwdef struct TermParams <: Params
	dims = (8, 8, 8)
	#β = [i/2 for i in 1.0:15.0]
	β = 1.0
	latticestart = :cold
	#sp2type = Sp2ElementB
	
	nterm = 50
	nnorm = 10 # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
	
	meanoffset = 5 # from which iteration to start calculating the mean

	nobs = 1 # measure observables every nobs cycles
end

Params(::Type{ObsParams}, x...) = ObsParams(x...) 
Params(::Type{TermParams}, x...) = TermParams(x...)

DrWatson.default_allowed(::TermParams) = (Real, String, Tuple)
DrWatson.allaccess(::TermParams) = (:dims, :β, :nterm, :nover, :nnorm)

