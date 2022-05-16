Base.@kwdef struct ObsParams <: Params
	save_plot    = false
	display_plot = true
	save_dat     = false
	save_jld2    = false
	save_df      = false
end

Base.@kwdef struct TermParams <: Params
	# dims = [(2, 8, 8, 8), (2, 10, 10, 10), (2, 12, 12, 12)]
	# β = sort(unique(vcat(6.43:0.01:6.49, 6.456:0.001:6.47)))	
	dims = [(2, 8, 8), (2, 10, 10), (2, 12, 12)]
	β = sort(unique(vcat(6.43:0.01:6.49, 6.456:0.001:6.47)))
	latticestart = :cold
	
	nterm = 300
	nnorm = 10 # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
	
	startobs = 100 # from which iteration to start taking measurements	
	nobs = 1 # measure observables every nobs cycles
	
	observables = (
		:avg_plaq	=> averageplaquette,
		:φ			=> expval_polyloop,
		:mod_φ		=> expval_modpolyloop,
		:χ			=> susceptibility,
		:χᵥ			=> susceptibility_pervolume,
		:S			=> action,
		:S²			=> actionsquared,
	)
end

DrWatson.default_allowed(::TermParams) = (Real, String, Tuple)
DrWatson.allaccess(::TermParams) = (:dims, :β, :nterm, :nover, :nnorm)

