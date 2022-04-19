Base.@kwdef struct ObsParams
	observables = ("avg_plaq" => averageplaquette, "action" => action)
	to_plot = true
	save_dat = true
	save_jld2 = true
	save_df = true
end

Base.@kwdef struct TermParams <: AbstractDict{String, Any}
	dims = [(8, 8, 8), (8, 8, 8, 8)]
	β = [i/2 for i in 1.0:15.0]
	latticestart = :cold
	sp2type = Sp2ElementB

	nterm = 100
	nnorm = 5 # after how many cycle to normalize lattice
	nover = 3 # how many cycles of overrelaxation to do
	
	meanoffset = 1 # from which iteration to start calculating the mean
end

function Base.iterate(X::TermParams, state)
	if length(state) > 0
		k = state[begin]
		return String(k) => getfield(X, k), Base.tail(state)
	end
	nothing
end
Base.iterate(X::TermParams) = Base.iterate(X, fieldnames(TermParams))
Base.length(::TermParams) = fieldcount(TermParams)
Base.get(X::TermParams, s::Symbol, ::Symbol) = Base.getproperty(X, s)
Base.get(X::TermParams, s, ::Symbol) = Base.getproperty(X, Symbol(s))

function DrWatson.dict_list(p::TermParams)
	[TermParams([d[String(field)] for field in fieldnames(TermParams)]...)
	for d in DrWatson.dict_list(Dict(p))]
end

DrWatson.default_allowed(::TermParams) = (Real, String, Tuple)
DrWatson.allaccess(::TermParams) = (:dims, :β, :nterm, :nover, :nnorm)