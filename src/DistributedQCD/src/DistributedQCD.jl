module DistributedQCD

using Distributed, DistributedArrays, StaticArrays, LinearAlgebra, Distributions
using Base: tail
using BlockArrays: mortar

export @everywhere, procs, workers, nworkers, initprocs, with_workers, @maybe_threaded
export SU2, asmatrix, normalizeSU2, normalizeSU2det
export Sp2, asmatrix, normalizeSp2
export Mask, Indices, Site, LocalLattice, Lattice, newlattice, evenmask, oddmask
export one_termalization!
export ObsConfig, averageplaquette, polyloop, corr_polyloop, expval_polyloop, expval_modpolyloop
export susceptibility, χ, susceptibility_pervolume, χᵥ, binder, gᵣ, action, actionsquared

"""
	initprocs(p; threads = 1, kwargs...)
Call `addprocs` with `p` processes passing the flags `--project --threads=...`, where the number of threads is decided by the `threads` argument.
"""
function initprocs(p; threads = 1, kwargs...) 
    ENV["JULIA_NODES"] = length(p)
    addprocs(p; exeflags = ["--project", "--threads=$threads"], kwargs...)
end

"""
	with_workers(f, args...; procs = workers())
Spawn a task calling `f(w, args...)` for each process in `workers()`. Useful with the `do end` statement to run a function on all processes exactly once.
`procs` is the list of processes that run `f`, defaults to all processes. The main worker will wait for all \
processes to finish before continuing.
"""
with_workers(f, args...; procs = workers()) = @sync for w in procs
	@spawnat w f(args...)
end

"""
This macro is identical to `Threads.@threads` except it doens't use threads if `Threads.nthreads() == 1`. \
This is used to avoid the overhead of multithreading in case only one thread is used.
"""
macro maybe_threaded(ex)
    if Threads.nthreads() == 1
        return esc(ex)
    else
        return esc(:(Threads.@threads $ex))
    end
end

include("SU2.jl")
include("Sp2.jl")
include("Lattice.jl")
include("CabibboMarinari.jl")
include("Observables.jl")

end # module DistributedQCD