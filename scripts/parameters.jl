abstract type Params <: AbstractDict{String, Any} end

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
	sp2type = nothing

	nterm = 50
	nnorm = 10 # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
	
	meanoffset = 5 # from which iteration to start calculating the mean

	nobs = 1 # measure observables every nobs cycles
end


#* ===== FUNCTIONS =====

Params(::Type{ObsParams}, x...) = ObsParams(x...) 
Params(::Type{TermParams}, x...) = TermParams(x...)

DrWatson.default_allowed(::TermParams) = (Real, String, Tuple)
DrWatson.allaccess(::TermParams) = (:dims, :β, :nterm, :nover, :nnorm)

function Base.iterate(X::Params, state)
	if length(state) > 0
		k = state[begin]
		return String(k) => getfield(X, k), Base.tail(state)
	end
	nothing
end
Base.iterate(X::Params) = Base.iterate(X, fieldnames(typeof(X)))
Base.length(X::Params) = fieldcount(typeof(X))
Base.get(X::Params, s::Symbol, ::Symbol) = Base.getproperty(X, s)
Base.get(X::Params, s, ::Symbol) = Base.getproperty(X, Symbol(s))

DrWatson.dict_list(p::Params) = [
	Params(
		typeof(p), 
		[d[String(field)] for field in fieldnames(typeof(p))]...
	)
	for d in DrWatson.dict_list(Dict(p))
]
