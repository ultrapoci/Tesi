import DrWatson, Statistics, DataFrames, Measurements, CSV, DelimitedFiles, ProgressMeter, Plots

function DrWatson._wsave(filename, data::Dict)
	if splitext(filename)[2] == ".dat" 
		DelimitedFiles.writedlm(filename, data, " = ")
	else
		DrWatson.save(filename, data)
	end
end

DrWatson._wsave(filename, data::DataFrames.DataFrame) = CSV.write(filename, data)
DrWatson._wsave(filename, p::Plots.Plot; kwargs...) = Plots.savefig(p, filename; kwargs...)

showall(x) = begin show(stdout, "text/plain", x); println() end

takeobservables(x::Pair) = (x.first,), (x.second,)
takeobservables(x) = first.(x), last.(x)

incremental_measurement(v) = [Measurements.measurement(Statistics.mean(v[1:i]), Statistics.std(v[1:i])) for i in 1:length(v)]

getpbar(n; desc = "Progress: ", enabled = true) = ProgressMeter.Progress(n, dt = 1, desc = desc, enabled = enabled, showspeed = true)
generate_showvalues(pairs...) = () -> [Tuple.(pairs)...]

totaliter(allparams) = sum((p[:nterm] for p in DrWatson.dict_list(allparams)))

Measurements.measurement(v::Vector) = Measurements.measurement(Statistics.mean(v), Statistics.std(v) / sqrt(length(v)))

function partition_workers(w, n; strategy = :atleast)
	par = convert.(Array, collect(Iterators.partition(w, n)))
	l = length(par[begin:end-1])
	if strategy == :atleast
		if length(par[end]) < length(par[begin])
			for (i, elem) in enumerate(par[end])
				push!(par[mod1(i, l)], elem)
			end
			par[begin:end-1]
		else
			par
		end
	elseif strategy == :atmost
		i = firstindex(par)
		while length(par[end]) < length(par[begin])
			push!(par[end], pop!(par[mod1(i, l)]))
			i += 1
		end
		par
	else
		throw(ArgumentError("strategy must be :atleast or :atmost, got :$strategy."))
	end
end

"""
	readinfo(f)
Read info from a CSV file `f`. Info is stored in the first two lines of the file (header + one line of data).
"""
function readinfo(f)
	d = CSV.read(f, DataFrames.DataFrame, limit = 1)[1, :]
	Dict(Symbol.(names(d)) .=> values(d))
end

"""
	readdata(f)
Read data from a CSV file `f`. Data is stored from line 3 (the header).
"""
readdata(f) = CSV.read(f, DataFrame, header = 3)

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

