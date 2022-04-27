function averageplaquette(L::Lattice{D}, inds::Indices{D}; log = false) where D
	log && @info "Measuring average plaquette..."
	current_workers = [w for w in procs(L)]
	n_current_workers = length(current_workers)
	# initialize a vector with nworkers() elements, of which each worker owns one cell
	partial = dzeros((n_current_workers,), current_workers, [n_current_workers])
	with_workers() do _	
		tot = 0.0
		for x in CartesianIndices(L[:L]), u in 1:D-1, v in u+1:D
			tot += tr(L[:L][x][u] * staple(L, v, u, inds[:L][x]))
		end
		partial[:L][begin] = tot
	end
	res = sum(partial)
	V = length(L) # lattice's volume
	np = D * (D - 1) ÷ 2 # number of plaquettes per site

	# the 4 term is due to 2N for N=2
	res / (4 * np * V) 
end
averageplaquette(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}; kwargs...) where D = averageplaquette(T[1], T[3]; kwargs...)
averageplaquette(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}; kwargs...) where D = averageplaquette(T.lattice, T.inds; kwargs...)

function all_ploops(L::Lattice{D}) where D
	current_workers = [procs(L)[i] for i in eachindex(procs(L))] # vector of workers that own L
	dist = size(procs(L)) # how are partitions distributed among workers

	# initialize an empty array of array: each cell belongs to a worker, and it contains a D-1 dimensional array to store partial product of time-like link 
	partial = DArray(I -> Array{Array{Sp2, D-1}}(undef, length.(I)...), dist, current_workers, dist);
	
	with_workers(procs = current_workers) do _
		localspacedims = tail(size(L[:L])) # size of the space part of local L
		# product along time direction for each space point
		time_direction = [l[1] for l in L[:L]] # links pointing in time direction
		partial[:L][begin] = [prod(time_direction[:, x]) for x in CartesianIndices(localspacedims)] # iterate through spatial indices of L[:L]
	end
	loc = convert(Array, partial) # bring partial into local process
	spacedist = tail(dist) # distrbutions of the space indices

	# product of each matrix of space points along the time axes. This is done because the time axes could be 
	# divide into different processes, and we need to multiply every partial polyakov loop.
	# w will be an array of arrays, one for each process, and each of these arrays
	# contains the product of time-link link in each lattice's position belonging to that process.
	# "mortar" transforms the array of arrays into a single array that corresponds to the spatial slice of L
	tr.(convert(Array, mortar([reduce(.*, loc[:, x]) for x in CartesianIndices(spacedist)])))
end
all_ploops(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}; kwargs...) where D = all_ploops(T[1]; kwargs...)
all_ploops(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}; kwargs...) where D = all_ploops(T.lattice; kwargs...)

function expval_ploop(L::Lattice{D}; log = false) where D
	log && @info "Measuring exp. value of Polyakov loops..."

	res = sum((all_ploops(L))) # sum of the traces of all polyakov loops
	Vₛ = prod(tail(size(L))) # volume of the space slice of the lattice L

	res / Vₛ
end
expval_ploop(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}; kwargs...) where D = expval_ploop(T[1]; kwargs...)
expval_ploop(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}; kwargs...) where D = expval_ploop(T.lattice; kwargs...)

function expval_modploop(L::Lattice{D}; log = false) where D
	log && @info "Measuring exp. value of Polyakov loops..."

	res = sum(abs.(all_ploops(L))) # sum of the absolute value of the traces of all polyakov loops
	Vₛ = prod(tail(size(L))) # volume of the space slice of the lattice L

	res / Vₛ
end
expval_modploop(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}; kwargs...) where D = expval_modploop(T[1]; kwargs...)
expval_modploop(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}; kwargs...) where D = expval_modploop(T.lattice; kwargs...)

function ploop(L::Lattice{D}, x::NTuple{Dm1, Int}; log = false) where {D, Dm1}
	Dm1 ≠ D-1 && throw(TypeError(ploop, NTuple{D-1, Int}, NTuple{Dm1, Int}))

	log && @info "Measuring Polyakov loop at $x..."

	tr(all_ploops(L)[x...])
end
ploop(L, x; kwargs...) = ploop(L, Tuple(x); kwargs...)
ploop(L, x...; kwargs...) = ploop(L, Tuple(x); kwargs...)
ploop(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, x; kwargs...) where D = ploop(T[1], x; kwargs...)
ploop(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, x; kwargs...) where D = ploop(T.lattice, x; kwargs...)

function corr_ploop(L::Lattice{D}, x::NTuple{Dm1, Int}, y::NTuple{Dm1, Int}; log = false) where {D, Dm1}
	Dm1 ≠ D-1 && throw(TypeError(ploop, NTuple{D-1, Int}, NTuple{Dm1, Int}))

	log && @info "Measuring two point correlation function of Polyakov loops at $x and $y..."

	loops = all_ploops(L)
	tr(loops[x...]) * tr(loops[y...])
end
corr_ploop(L, x, y; kwargs...) = corr_ploop(L, Tuple(x), Tuple(y); kwargs...)
corr_ploop(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, x, y; kwargs...) where D = corr_ploop(T[1], x, y; kwargs...)
corr_ploop(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, x, y; kwargs...) where D = corr_ploop(T.lattice, x, y; kwargs...)

χ = susceptibility
function susceptibility(L::Lattice{D}; log = false) where D
	log && @info "Measuring susceptibility..."

	sum(all_ploops(L) .^ 2)
end
susceptibility(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}; kwargs...) where D = susceptibility(T[1]; kwargs...)
susceptibility(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}; kwargs...) where D = susceptibility(T.lattice; kwargs...)

gᵣ = binder
function binder(L::Lattice{D}; log = false) where D
	log && @info "Measuring Binder cumulant..."

	loops = all_ploops(L)
	Vₛ = length(loops) # spatial volume
	φ⁴ = sum(loops .^ 4)
	φ² = sum(loops .^ 2)

	Vₛ * φ⁴ / (φ²^2) - 3
end
binder(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}; kwargs...) where D = binder(T[1]; kwargs...)
binder(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}; kwargs...) where D = binder(T.lattice; kwargs...)

#=
function ploop(L::Lattice{D}, x::NTuple{Dm1, Int}; log = false) where {D, Dm1}
	Dm1 ≠ D-1 && throw(TypeError(corr_loop, CartesianIndex{D-1}, CartesianIndex{Dm1}))

	log && @info "Measuring Polyakov loop at $x..."

	current_workers = [procs(L)[i] for i in eachindex(procs(L))] # vector of workers that own L
	dist = size(procs(L)) # how are partitions distributed among workers
	X = distribute(Array{Union{Sp2, Missing}}(undef, dist...), procs = current_workers, dist = [dist...])
	
	with_workers(procs = current_workers) do _
		if all(x .∈ tail(localindices(L)))
			time_direction = [l[1] for l in L[:L]]
			X[:L][begin] = prod(time_direction[:, x...])
		end
	end
	loc = convert(Array, X)
	spacedist = tail(dist)

	# only one spatial index will have a time slice of Sp2 matrices that are not missing
	for x in CartesianIndices(spacedist)
		if !ismissing(loc[1, x])
			return tr(reduce(*, loc[:, x]))
		end
	end
end
ploop(L, x::CartesianIndex; kwargs...) = ploop(L, Tuple(x); kwargs...)
ploop(L, x::Vararg{Int}; kwargs...) = ploop(L, Tuple(x); kwargs...)
ploop(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, x; kwargs...) where D = ploop(T[1], x; kwargs...)
ploop(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, x; kwargs...) where D = ploop(T.lattice, x; kwargs...)
=#