Base.@kwdef struct ObsParams <: Params
	save_plot    = false
	display_plot = true
	save_dat     = false
	save_jld2    = false
	save_df      = false
end

Base.@kwdef struct TermParams <: Params
	β = 7.635
	dims = (20, 16, 16, 16)
	latticestart = :cold
	
	nterm = 1_000
	nnorm = 100 # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
	
	observables = (
		:plaqs	=> plaquettesum,
		:φ		=> polyloop_sum,
		:φ²		=> polyloop2_sum,
		:φ⁴		=> polyloop4_sum,
		:modφ	=> polyloopmod_sum,
		:R2		=> (C::ObsConfig; kwargs...) -> corr_polyloop(2, C; kwargs...),
		:R3		=> (C::ObsConfig; kwargs...) -> corr_polyloop(3, C; kwargs...),
		:R4		=> (C::ObsConfig; kwargs...) -> corr_polyloop(4, C; kwargs...),
		:R5		=> (C::ObsConfig; kwargs...) -> corr_polyloop(5, C; kwargs...),
		:R6		=> (C::ObsConfig; kwargs...) -> corr_polyloop(6, C; kwargs...),
		:R7		=> (C::ObsConfig; kwargs...) -> corr_polyloop(7, C; kwargs...),
		:R8		=> (C::ObsConfig; kwargs...) -> corr_polyloop(8, C; kwargs...),
		# :V		=> (C::ObsConfig; kwargs...) -> length(C.L),
		# :Nₜ		=> (C::ObsConfig; kwargs...) -> first(size((C.L))),
		# :Vₛ		=> (C::ObsConfig; kwargs...) -> spatialvolume(C.L),
	)
	
	# startobs = 400 # from which iteration to start taking measurements	
	# nobs = 1 # measure observables every nobs cycles
	
	#= observables = (
		:avg_plaq	=> averageplaquette,
		:φ			=> expval_polyloop,
		:mod_φ		=> expval_modpolyloop,
		:φ²			=> expval_polyloop_squared,
		:χ			=> susceptibility,
		:χᵥ			=> susceptibility_pervolume,
		:S			=> action,
		:S²			=> actionsquared,
		:gᵣ			=> binder
	) =#
end

DrWatson.default_allowed(::TermParams) = (Real, String, Tuple)
DrWatson.allaccess(::TermParams) = (:dims, :β, :nterm, :nover, :nnorm)

