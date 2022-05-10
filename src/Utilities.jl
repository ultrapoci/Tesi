import DrWatson, Statistics, DataFrames, Measurements, CSV, DelimitedFiles, ProgressMeter

function DrWatson._wsave(filename, data::Dict)
	if splitext(filename)[2] == ".dat" 
		DelimitedFiles.writedlm(filename, data, " = ")
	else
		DrWatson.save(filename, data)
	end
end

DrWatson._wsave(filename, data::DataFrames.DataFrame) = CSV.write(filename, data)

showall(x) = begin show(stdout, "text/plain", x); println() end

takeobservables(x::Pair) = (x.first,), (x.second,)
takeobservables(x) = first.(x), last.(x)

keys_to_symbol(d::Dict) = Dict(Symbol.(keys(d)) .=> values(d))
keys_to_string(d::Dict) = Dict(String.(keys(d)) .=> values(d))

incremental_measurement(v) = [Measurements.measurement(Statistics.mean(v[1:i]), Statistics.std(v[1:i])) for i in 1:length(v)]

getpbar(n; desc = "Progress: ", enabled = true) = ProgressMeter.Progress(n, dt = 1, desc = desc, enabled = enabled, showspeed = true)
generate_showvalues(pairs...) = () -> [Tuple.(pairs)...]

totaliter(allparams) = sum((p[:nterm] for p in DrWatson.dict_list(allparams)))


#* ===== PARAMETERS =====

abstract type Params <: AbstractDict{Symbol, Any} end

function Base.iterate(X::Params, state)
	if length(state) > 0
		k = state[begin]
		return k => getfield(X, k), Base.tail(state)
	end
	nothing
end
Base.iterate(X::T) where T <: Params = Base.iterate(X, fieldnames(T))
Base.length(::T) where T <: Params = fieldcount(T)
Base.get(X::Params, s::Symbol, ::Symbol) = Base.getproperty(X, s)
Base.get(X::Params, s, ::Symbol) = Base.getproperty(X, Symbol(s))

DrWatson.dict_list(p::T) where T <: Params = [T(; d...)	for d in DrWatson.dict_list(Dict(p))]

