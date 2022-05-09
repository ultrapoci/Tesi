Base.@kwdef struct ObsParams <: Params
	save_plot    = false
	display_plot = true
	save_dat     = false
	save_jld2    = false
	save_df      = false
end

Base.@kwdef struct TermParams <: Params
	dims = (2, 8, 8, 8)
	#β = [i for i in 1.0:15.0]
	β = 6.45
	latticestart = :cold
	sp2type = try Sp2ElementB catch; missing end
	
	nterm = 400
	nnorm = 10 # after how many cycle to normalize lattic40, 40e
	nover = 3 # how many cycles of overrelaxation to do
	
	startobs = 100 # from which iteration to start taking measurements	
	nobs = 1 # measure observables every nobs cycles
	
	observables = (
		:avg_plaq		=> averageplaquette,
		:polyloop		=> expval_polyloop,
		:mod_polyloop	=> expval_modpolyloop,
		:χ_product		=> susceptibility,
		:χᵥ				=> susceptibility_pervolume,
		:action			=> action,
		:action_squared	=> actionsquared,
		#:χ_product		=> susceptibility2,
		#:χᵥ_product	=> susceptibility_pervolume2,
	)
end

DrWatson.default_allowed(::TermParams) = (Real, String, Tuple)
DrWatson.allaccess(::TermParams) = (:dims, :β, :nterm, :nover, :nnorm)

