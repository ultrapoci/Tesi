using DrWatson, Statistics, ProgressMeter, DataFrames
import CSV, DelimitedFiles

function DrWatson._wsave(filename, data::Dict)
	if splitext(filename)[2] == ".dat" 
		DelimitedFiles.writedlm(filename, data, " = ")
	else
		save(filename, data)
	end
end

DrWatson._wsave(filename, data::DataFrame) = CSV.write(filename, data)

showall(x) = begin show(stdout, "text/plain", x); println() end

takeobservables(x::Pair) = (x.first,), (x.second,)
takeobservables(x) = first.(x), last.(x)

to_symbol(d) = Dict(Symbol.(keys(d)) .=> values(d))

function incrementalmean(v, offset::Integer = 1)
	if offset ∉ 1:length(v)
		throw(ArgumentError("Given offset = $offset must be positive and less than or equal to length(v) = $(length(v))"))
	end

	[mean(v[offset:i]) for i in offset:length(v)], [std(v[offset:i]) for i in offset:length(v)]
end

getpbar(n; desc = "Progress: ", enabled = true) = Progress(n, dt = 1, desc = desc, enabled = enabled, showspeed = true)
generate_showvalues(pairs...) = () -> [Tuple.(pairs)...]


#* ===== PARAMETERS =====

abstract type Params <: AbstractDict{String, Any} end

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

