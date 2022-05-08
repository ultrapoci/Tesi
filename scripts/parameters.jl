Base.@kwdef struct ObsParams <: Params
	save_plot    = false
	display_plot = true
	save_dat     = false
	save_jld2    = false
	save_df      = false
end

Base.@kwdef struct TermParams <: Params
	dims = (2, 40, 40)
	#β = [i for i in 1.0:15.0]
	β = [10.3, 10.4, 10.5]
	latticestart = :cold
	sp2type = try Sp2ElementB catch; missing end
	
	nterm = 15_000
	nnorm = 100 # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
	
	startobs = 300 # from which iteration to start taking measurements	
	nobs = 1 # measure observables every nobs cycles
	
	observables = (
		:avg_plaq			=> averageplaquette,
		:polyloop			=> expval_polyloop,
		:mod_polyloop		=> expval_modpolyloop,
		#"χ_sum_squared"  => susceptibility,
		#"χ_product"      => susceptibility2,
		#"χᵥ_sum_squared" => susceptibility_pervolume,
		#"χᵥ_product" 	 => susceptibility_pervolume2,
		#"action"         => action,
		#"action_squared" => actionsquared,
	)
end

DrWatson.default_allowed(::TermParams) = (Real, String, Tuple)
DrWatson.allaccess(::TermParams) = (:dims, :β, :nterm, :nover, :nnorm)

