using DrWatson
@quickactivate "Tesi"

if "TERMALIZATION_LOG" ∉ keys(ENV)
	ENV["TERMALIZATION_LOG"] = "0"
end

if "TERMALIZATION_NPROCS" ∉ keys(ENV)
	ENV["TERMALIZATION_NPROCS"] = "4"
end

using DistributedQCD, DataFrames, Suppressor, ProgressMeter
nworkers() == 1 && initprocs(parse(Int, ENV["TERMALIZATION_NPROCS"]))
@everywhere using DistributedQCD

include(srcdir("Utilities.jl"))
try 
	includet(scriptsdir("parameters.jl")) 
catch 
	@warn "Including parameters without Revise.jl"
	include(scriptsdir("parameters.jl")) 
end


#* ===== TERMALIZATION =====

function termalization!(L, params; log = false)
	@unpack β, nterm, nover, nnorm = params

	pbar = getpbar(nterm, desc = " Distributed Termalization", enabled = !log)

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		next!(pbar, showvalues = generate_showvalues(:iter => n, :total => nterm))
	end
end

function termalization(params; procs = workers(), kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart, procs = procs)
	termalization!(L, params; kwargs...)
	L
end

function termalization!(L, params, observable::Function, v::Vector; log = false)
	@unpack β, nterm, nover, nnorm, nobs, startobs = params

	pbar = getpbar(nterm, desc = " Distributed Termalization", enabled = !log)

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		T = (L..., β = β) # used to pass β to observables that require it
		if n ≥ startobs && n % nobs == 0
			push!(v, observable(T; log = log, iter = n))
		end
		next!(pbar, showvalues = generate_showvalues(:iter => n, :total => nterm))
	end
	
	# make sure the last iteration is measured
	nterm % nobs ≠ 0 && push!(v, observable(L; log = log, iter = n))
end

function termalization(params, observable::Function; procs = workers(), kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart, procs = procs)
	v = Float64[]
	termalization!(L, params, observable, v; kwargs...)
	v, L
end

function termalization!(L, params, observables, v; log = false)
	@unpack β, nterm, nover, nnorm, startobs, nobs = params

	pbar = getpbar(nterm, desc = " Distributed Termalization", enabled = !log)

	for n in 1:nterm
		log && @info "Termalization" n nterm
		one_termalization!(L, nover, β, n % nnorm == 0; log = log, iter = n)
		T = (L..., β = β) # used to pass β to observables that require it
		if n ≥ startobs && n % nobs == 0
			push!(v, [obs(T; log = log, iter = n) for obs in observables])
		end
		next!(pbar, showvalues = generate_showvalues(:iter => n, :total => nterm))
	end

	# make sure the last iteration is measured
	nterm % nobs ≠ 0 && push!(v, [obs(L; log = log, iter = n) for obs in observables])
end

function termalization(params, observables; procs = workers(), kwargs...)
	@unpack dims, latticestart = params
	L = newlattice(dims..., start = latticestart, procs = procs)
	v = Vector{Float64}[]
	termalization!(L, params, observables, v; kwargs...)
	reduce(hcat, v)', L
end


#* ===== RUN =====

function run(allparams, folder = ""; save = true)
	!save && @warn "Current simulation is not going to be saved!"

	@unpack observables = allparams
	obsnames, obsfunctions = takeobservables(observables)
	dicts = Dict[]

	for params in dict_list(allparams)
		@unpack startobs = params
		display(params)
		d = Dict(params)
		obsmeasurements, L = termalization(params, obsfunctions, log = ENV["TERMALIZATION_LOG"] == "1")
		d[:data] = DataFrame(obsmeasurements, [obsnames...])
		d[:L] = (lattice = convert(Array, L.lattice), mask = convert(Array, L.mask), inds = convert(Array, L.inds))
		if save
			jld2name = savename(params, "jld2", sort = false)
			@info "Saving jld2 file '$jld2name' in $(datadir(folder))..."
			@suppress_err safesave(datadir(folder, jld2name), d) # suppress warnings about symbols converted to strings
		end
		push!(dicts, d)
		println("="^displaysize(stdout)[2])
	end

	dicts
end

run(; kwargs...) = run(TermParams(); kwargs...)
run(s::String; kwargs...) = run(TermParams(), s; kwargs...)