Base.@kwdef struct ObsParams <: Params
	save_plot    = false
	display_plot = true
	save_dat     = false
	save_jld2    = false
	save_df      = false
end

Base.@kwdef struct TermParams <: Params
	dims = [(4, L, L) for L in 20:10:80]
	β = [@onlyif(:dims==(4, 80, 80), [6.5, 8.0, 8.5])...,
		sort(unique(vcat(21.0:1.0:30.0)))...]
	# β = sort(unique(vcat(
	# 	4.0:1.0:20.0,
	# 	6.0:0.5:11.0
	# )))
	latticestart = :cold
	
	nterm = 15_000
	nnorm = 100 # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
	
	startobs = 400 # from which iteration to start taking measurements	
	nobs = 1 # measure observables every nobs cycles
	
	observables = (
		:avg_plaq	=> averageplaquette,
		:φ			=> expval_polyloop,
		:mod_φ		=> expval_modpolyloop,
		:φ²			=> expval_polyloop_squared,
		:χ			=> susceptibility,
		:χᵥ			=> susceptibility_pervolume,
		:S			=> action,
		:S²			=> actionsquared,
		:gᵣ			=> binder
	)
end

DrWatson.default_allowed(::TermParams) = (Real, String, Tuple)
DrWatson.allaccess(::TermParams) = (:dims, :β, :nterm, :nover, :nnorm)

