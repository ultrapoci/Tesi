using DrWatson, Statistics, ProgressMeter, DataFrames, Measurements
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

incremental_measurement(v) = [measurement(mean(v[1:i]), std(v[1:i])) for i in 1:length(v)]

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
Base.iterate(X::T) where T <: Params = Base.iterate(X, fieldnames(T))
Base.length(::T) where T <: Params = fieldcount(T)
Base.get(X::Params, s::Symbol, ::Symbol) = Base.getproperty(X, s)
Base.get(X::Params, s, ::Symbol) = Base.getproperty(X, Symbol(s))

DrWatson.dict_list(p::T) where T <: Params = [
	T([d[String(field)] for field in fieldnames(T)]...)
	for d in DrWatson.dict_list(Dict(p))
]

