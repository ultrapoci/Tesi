module DistributedQCD

using Distributed, DistributedArrays, StaticArrays, LinearAlgebra, Distributions

export @everywhere, workers, nworkers, initprocs, with_workers

"""
	initprocs(n; kwargs...)
Calls `addprocs` with `n` processes and the flag `exeflags="--project"`. All other `kwargs` are passed to `addprocs`. 
"""
initprocs(n; kwargs...) = addprocs(n, exeflags="--project", kwargs...)

"""
	with_workers(f)
Spawn a task calling `f()` for each process in `workers()`. Useful with the `do end` statement to run a function on all processes exactly once.
`f` takes the id of the worker as the only parameter.
"""
with_workers(f) = @sync for w in workers()
	@spawnat w f(w)
end

include("SU2.jl")
include("Sp2.jl")
include("Lattice.jl")
include("CabibboMarinari.jl")
include("Observables.jl")

end # module StaticQCD