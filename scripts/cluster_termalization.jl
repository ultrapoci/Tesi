using DrWatson
@quickactivate "Tesi"

if "TERMALIZATION_LOG" ∉ keys(ENV)
	ENV["TERMALIZATION_LOG"] = "0"
end

using Distributed, DistributedQCD, DataFrames, Suppressor, ProgressMeter

include(srcdir("Utilities.jl"))
try 
	using Revise
	includet(scriptsdir("parameters.jl")) 
catch 
	@warn "Including parameters without Revise.jl"
	include(scriptsdir("parameters.jl")) 
end


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
	@unpack β, nterm, nover, nnorm, nobs, startobs = params

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		if n ≥ startobs && n % nobs == 0
			push!(v, observable(ObsConfig(L, β); log = log, iter = n))
		end
		put!(channel, true)
	end
	
	# make sure the last iteration is measured
	nterm % nobs ≠ 0 && push!(v, observable(ObsConfig(L, β); log = log, iter = n))
end

function termalization(params, observable::Function, channel; procs = workers(), kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart, procs = procs)
	v = Float64[]
	termalization!(L, params, observable, v, channel; kwargs...)
	v, L
end

function termalization!(L, params, observables, v, channel; log = false)
	@unpack β, nterm, nover, nnorm, startobs, nobs = params

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		if n ≥ startobs && n % nobs == 0
			push!(v, [obs(ObsConfig(L, β); log = log, iter = n) for obs in observables])
		end
		put!(channel, true)
	end

	# make sure the last iteration is measured
	nterm % nobs ≠ 0 && push!(v, [obs(ObsConfig(L, β); log = log, iter = n) for obs in observables])
end

function termalization(params, observables, channel; procs = workers(), kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart, procs = procs)
	v = Vector{Float64}[]
	termalization!(L, params, observables, v, channel; kwargs...)
	reduce(hcat, v)', L
end


#* ===== RUN =====

function run(allparams; n = nothing, strategy = :atleast, folder = "", save = true)
	strategy ∉ [:atleast, :atmost] && throw(ArgumentError("strategy must be :atleast or :atmost, got :$strategy."))
	!save && @warn "Current simulation is not going to be saved!"

	# setup workers pool
	wp, pool, pbardesc = if isnothing(n)
		@info "Dividing workers according to the node they belong."
		nodes = parse(Int, ENV["JULIA_NODES"])
		p = nodes == 1 ? [workers()] : procs.(workers()[begin:nodes])
		WorkerPool(workers()[begin:nodes]), p, "Termalization using $nodes nodes..."
	else
		p = partition_workers(workers(), n, strategy = strategy)
		WorkerPool(first.(p)), p, "Termalization using $(length(p)) groups of workers..."
	end
	getpool(id) = pool[findfirst(x -> id ∈ x, pool)]

	# setup progress bar
	# pbar = ProgressBar(expand = true, columns = :detailed, refresh_rate = 1)
	# job = addjob!(pbar, N = totaliter(allparams), description = pbardesc)
	total_iter = totaliter(allparams)
	pbar = getpbar(total_iter, desc = pbardesc, enabled = ENV["TERMALIZATION_LOG"] != "1")
	channel = RemoteChannel(()->Channel{Bool}(), 1)
	@async begin
		i = 0 
		while take!(channel)
			i += 1
			next!(pbar, showvalues = generate_showvalues(:iter => i, :total => total_iter))
		end
	end
	
	logchannel = RemoteChannel(()->Channel{Tuple}(), 1)
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
	dicts = pmap(wp, dict_list(allparams)) do params
		put!(logchannel, (true, "From worker $(myid()): starting new termalization with workers = $(getpool(myid()))"))
		d = Dict(params)
		obsmeasurements, L = termalization(params, obsfunctions, channel, procs = getpool(myid()), log = ENV["TERMALIZATION_LOG"] == "1")
		d[:data] = DataFrame(obsmeasurements, collect(obsnames))
		d[:L] = (lattice = convert(Array, L.lattice), mask = convert(Array, L.mask), inds = convert(Array, L.inds))
		if save
			jld2name = savename(params, "jld2", sort = false)
			put!(logchannel, (true, "From worker $(myid()): saving jld2 file '$jld2name' in $(datadir(folder))"))
			@suppress_err safesave(datadir(folder, jld2name), d) # suppress warnings about symbols converted to strings
		end
		d
	end
	put!(channel, false)
	put!(logchannel, (false,))
	dicts
end
run(; kwargs...) = run(TermParams(); kwargs...)

#=
function oldrun(allparams, folder = ""; save = true)
	!save && @warn "Current simulation is not going to be saved!"

	nodes = parse(Int, ENV["JULIA_NODES"])
	wp = WorkerPool(workers()[begin:nodes])
	pbar = Progress(totaliter(allparams), dt = 1)
	channel = RemoteChannel(()->Channel{Bool}(), 1)

	@async while take!(channel)
		next!(pbar)
	end

	@unpack observables = allparams
	obsnames, obsfunctions = takeobservables(observables)

	dicts = pmap(wp, dict_list(allparams)) do params
		@unpack startobs = params
		d = Dict(params)
		obsmeasurements, L = termalization(params, obsfunctions, channel, procs = procs(myid()), log = ENV["TERMALIZATION_LOG"] == "1")
		d[:data] = DataFrame(obsmeasurements, [obsnames...])
		d[:L] = (lattice = convert(Array, L.lattice), mask = convert(Array, L.mask), inds = convert(Array, L.inds))
		if save
			jld2name = savename(params, "jld2", sort = false)
			@info "Saving jld2 file $jld2name in $(datadir(folder))..."
			@suppress_err safesave(datadir(folder, jld2name), d) # suppress warnings about symbols converted to strings
		end
		d
	end
	put!(channel, false)
	dicts
end
=#