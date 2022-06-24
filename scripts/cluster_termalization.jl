using DrWatson
@quickactivate "Tesi"

if "TERMALIZATION_LOG" ∉ keys(ENV)
	ENV["TERMALIZATION_LOG"] = "0"
end

@info "Importing packages..."
using Distributed, DistributedQCD, DataFrames, Suppressor, ProgressMeter
import Base: tail

@info "Including 'Utilities.jl'..."
include(srcdir("Utilities.jl"))
try
	using Revise
	@info "Including 'parameters.jl'..."
	includet(scriptsdir("parameters.jl")) 
catch 
	@warn "Including parameters without Revise.jl"
	include(scriptsdir("parameters.jl")) 
end

@info "Building termalization functions..."
#* ===== TERMALIZATION =====

function termalization!(L, params, channel; log = false)
	@unpack β, nterm, nover, nnorm = params

	pbar = getpbar(nterm, desc = " Distributed Termalization", enabled = !log)

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		put!(channel, true)
	end
end

function termalization(params, channel; procs = workers(), kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart, procs = procs)
	termalization!(L, params, channel; kwargs...)
	L
end

function termalization!(L, params, observable::Function, v::Vector, channel; log = false)
	@unpack β, nterm, nover, nnorm = params

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		C = ObsConfig(L, β)
		push!(v, observable(C; log = log, iter = n))
		put!(channel, true)
	end
end

function termalization(params, observable::Function, channel; procs = workers(), kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart, procs = procs)
	v = Float64[]
	termalization!(L, params, observable, v, channel; kwargs...)
	v, L
end

function termalization!(L, params, observables, v, channel; log = false)
	@unpack β, nterm, nover, nnorm = params

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		C = ObsConfig(L, β)
		push!(v, [obs(C; log = log, iter = n) for obs in observables])
		put!(channel, true)
	end
end

function termalization(params, observables, channel; procs = workers(), kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart, procs = procs)
	v = Vector{Float64}[]
	termalization!(L, params, observables, v, channel; kwargs...)
	reduce(hcat, v)', L
end


#* ===== RUN =====

function run(allparams, n::Int; strategy = :atleast, folder = "", save = true, digits = 3)
	strategy ∉ (:atleast, :atmost) && throw(ArgumentError("strategy must be :atleast or :atmost, got :$strategy."))
	if save
		@info "Data will be saved in $(datadir(folder))" 
	else
		@warn "Current simulation is not going to be saved!"
	end
	
	#= wp, pool, pbardesc = if isnothing(n)
		@info "Dividing workers according to the node they belong to."
		nodes = parse(Int, ENV["JULIA_NODES"])
		p = nodes == 1 ? [workers()] : procs.(workers()[begin:nodes])
		WorkerPool(workers()[begin:nodes]), p, "$nodes nodes"
	else
		p = partition_workers(workers(), n, strategy = strategy)
		WorkerPool(first.(p)), p, "$(length(p)) groups of workers"
	end =#

	# setup workers pool
	pool = partition_workers(workers(), n, strategy = strategy)
	wp = WorkerPool(first.(pool))
	pbardesc = "$(length(pool)) groups of workers"
	getpool(id) = pool[findfirst(x -> id ∈ x, pool)]

	# setup progress bar
	total_iter = totaliter(allparams)
	pbar = getpbar(total_iter, desc = pbardesc, enabled = ENV["TERMALIZATION_LOG"] != "1")
	pbarchannel = RemoteChannel(()->Channel{Bool}(), 1)
	logchannel = RemoteChannel(()->Channel{Tuple}(), 1)
	@async begin
		i = 0 
		while take!(pbarchannel)
			i += 1
			next!(pbar, showvalues = generate_showvalues(:iter => i, :total => total_iter))
		end
	end	
	@async while true
		(b, s) = take!(logchannel)
		if b
			println(rpad("\r" * s, displaysize(stdout)[2]))
		else
			break
		end
	end

	@unpack observables = allparams
	obsnames, obsfunctions = takeobservables(observables)

	dicts = pmap(wp, sort(dict_list(allparams), by=(x->x.dims))) do params
		put!(logchannel, (true, rpad("Worker $(myid()):", 15) * "dims=$(params.dims), β=$(params.β) - with workers = $(getpool(myid()))"))
		
		d = Dict(params)
		@unpack dims = params
		d[:V] = prod(dims) # volume
		d[:Vₛ] = prod(tail(dims)) # spatial volume
		d[:Nₜ] = first(dims) # time length

		obsmeasurements, L = termalization(
			params, 
			obsfunctions, 
			pbarchannel, 
			procs = getpool(myid()), 
			log = ENV["TERMALIZATION_LOG"] == "1"
		)
		
		d[:data] = DataFrame(obsmeasurements, collect(obsnames))
		d[:L] = NamedTuple(keys(L) .=> convert.(Array, values(L))) # save L to dataframe
		
		if save
			jld2name = savename(params, "jld2", digits = digits, sort = false)
			put!(logchannel, (true, rpad("Worker $(myid()):", 15) * "saving '$jld2name' in $(datadir(folder))"))
			@suppress_err safesave(datadir(folder, jld2name), d) # suppress warnings about symbols converted to strings
		end
		
		d
	end

	put!(pbarchannel, false)
	put!(logchannel, (false,))

	dicts
end
run(n::Int; kwargs...) = run(TermParams(), n; kwargs...)