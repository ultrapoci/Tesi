Base.@kwdef struct ObsParams <: Params
	save_plot = false
	display_plot = true
	save_dat = false
	save_jld2 = false
	save_df = false
	observables = (
		#"avg_plaq" => averageplaquette,
		#"mod_polyloop" => expval_modpolyloop,
		"sum * V" => (L; kwargs...) -> sum(DistributedQCD.all_polyloops(L) .^ 2) * 512,
		"sum" => (L; kwargs...) -> sum(DistributedQCD.all_polyloops(L) .^ 2),
		"sum / V" => (L; kwargs...) -> sum(DistributedQCD.all_polyloops(L) .^ 2) / 512,
		"sum / V²" => (L; kwargs...) -> sum(DistributedQCD.all_polyloops(L) .^ 2) / 512^2,
	)
end

Base.@kwdef struct TermParams <: Params
	dims = (2, 8, 8, 8)
	#β = [i/2 for i in 1.0:15.0]
	β = 6.45 / 2
	latticestart = :cold
	sp2type = try Sp2ElementB catch; missing end
	
	nterm = 100
	nnorm = 10 # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
	
	meanoffset = 10 # from which iteration to start calculating the mean

	nobs = 1 # measure observables every nobs cycles
end

Params(::Type{ObsParams}, x...) = ObsParams(x...) 
Params(::Type{TermParams}, x...) = TermParams(x...)

DrWatson.default_allowed(::TermParams) = (Real, String, Tuple)
DrWatson.allaccess(::TermParams) = (:dims, :β, :nterm, :nover, :nnorm)

