Base.@kwdef struct ObsParams <: Params
	save_plot    = false
	display_plot = true
	save_dat     = false
	save_jld2    = false
	save_df      = false
end

Base.@kwdef struct TermParams <: Params
	dims = (8, 8, 8)
	β = [1.0, 2.0, 6.45]
	latticestart = :cold
	
	nterm = 20
	nnorm = 10 # after how many cycle to normalize lattic40, 40e
	nover = 3 # how many cycles of overrelaxation to do
	
	startobs = 10 # from which iteration to start taking measurements	
	nobs = 1 # measure observables every nobs cycles
	
	observables = (
		:avg_plaq		=> averageplaquette,
		:polyloop		=> expval_polyloop,
		:mod_polyloop	=> expval_modpolyloop,
		:χ				=> susceptibility,
		:χᵥ				=> susceptibility_pervolume,
		:action			=> action,
		:action_squared	=> actionsquared,
	)
end

DrWatson.default_allowed(::TermParams) = (Real, String, Tuple)
DrWatson.allaccess(::TermParams) = (:dims, :β, :nterm, :nover, :nnorm)

