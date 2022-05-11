using DrWatson
@quickactivate "Tesi"

if "TERMALIZATION_LOG" ∉ keys(ENV)
	ENV["TERMALIZATION_LOG"] = "0"
end

try
	using OhMyREPL, Revise
catch e
	@warn "Failed importing OhMyREPL, Revise: " e
end
using Distributed, DistributedQCD, DataFrames, Suppressor, Term.progress

include(srcdir("Utilities.jl"))
try 
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
		T = (L..., β = β) # used to pass β to observables that require it
		if n ≥ startobs && n % nobs == 0
			push!(v, observable(T; log = log, iter = n))
		end
		put!(channel, true)
	end
	
	# make sure the last iteration is measured
	nterm % nobs ≠ 0 && push!(v, observable(L; log = log, iter = n))
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
		T = (L..., β = β) # used to pass β to observables that require it
		if n ≥ startobs && n % nobs == 0
			push!(v, [obs(T; log = log, iter = n) for obs in observables])
		end
		put!(channel, true)
	end

	# make sure the last iteration is measured
	nterm % nobs ≠ 0 && push!(v, [obs(L; log = log, iter = n) for obs in observables])
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
		WorkerPool(workers()[begin:nodes]), procs.(workers()[begin:nodes]), "Termalization using $nodes nodes..."
	else
		p = partition_workers(workers(), n, strategy = strategy)
		WorkerPool(first.(p)), p, "Termalization using $(length(p)) groups of workers..."
	end
	getpool(id) = pool[id .∈ pool][begin]

	# setup progress bar
	pbar = ProgressBar(expand = true, columns = :detailed, refresh_rate = 1)
	job = addjob!(pbar, N = totaliter(allparams), description = pbardesc)
	channel = RemoteChannel(()->Channel{Bool}(), 1)
	@async while take!(channel)
		update!(job)
	end

	@unpack observables = allparams
	obsnames, obsfunctions = takeobservables(observables)
	start!(pbar)
	dicts = pmap(wp, dict_list(allparams)) do params
		@info "Starting new termalization: workers = $(getpool(myid()))"
		d = Dict(params)
		obsmeasurements, L = termalization(params, obsfunctions, channel, procs = getpool(myid()), log = ENV["TERMALIZATION_LOG"] == "1")
		d[:data] = DataFrame(obsmeasurements, collect(obsnames))
		d[:L] = (lattice = convert(Array, L.lattice), mask = convert(Array, L.mask), inds = convert(Array, L.inds))
		if save
			jld2name = savename(params, "jld2", sort = false)
			@info "Saving jld2 file '$jld2name' in $(datadir(folder))..."
			@suppress_err safesave(datadir(folder, jld2name), d) # suppress warnings about symbols converted to strings
		end
		d
	end
	put!(channel, false)
	stop!(pbar)
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