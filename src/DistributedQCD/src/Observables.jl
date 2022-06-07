struct ObsConfig{D}
	L::Lattice{D}
	inds::Indices{D}
	β::Real
	polyloops::Array{Float64}
	plaquette_sum::Real

	ObsConfig(L::Lattice{D}, inds::Indices{D}, β::Real) where D = new{D}(L, inds, β, all_polyloops(L), plaquettesum(L, inds))
	ObsConfig(T::NamedTuple, β::Real) = new{ndims(T.lattice)}(T.lattice, T.inds, β, all_polyloops(T.lattice), plaquettesum(T.lattice, T.inds))
end

#* ===== average plaquette =====
"""
	averageplaquette(L::Lattice{D}, inds::Indices{D}; log = false, iter = missing) where D
	averageplaquette(T::NamedTuple; kwargs...)
	averageplaquette(C::ObsConfig; kwargs...)
Return the value of the average plaquette of the lattice `L` divided by ``2N``, where ``N`` is the dimensions of ``Sp(N)``. In this case, ``N = 2``.
"""
averageplaquette(L::Lattice{D}, inds::Indices{D}; kwargs...) where D = _averageplaquette(plaquettesum(L, inds), length(L), D; kwargs...)
averageplaquette(T::NamedTuple; kwargs...) = averageplaquette(T.lattice, T.inds; kwargs...)
averageplaquette(C::ObsConfig{D}; kwargs...) where D = _averageplaquette(C.plaquette_sum, length(C.L), D; kwargs...)

function _averageplaquette(plaquette_sum::Real, volume::Int, D::Int; log = false, iter = missing)
	log && @info "Measuring average plaquette..." iter
	np = D * (D - 1) ÷ 2 # number of plaquettes per site
	plaquette_sum / (4 * np * volume) # the 4 term is due to 2N for N=2
end


#* ===== action =====
"""
	action(L::Lattice{D}, inds::Indices{D}; β, log = false, iter = missing) where D
	action(T::NamedTuple; kwargs...)
	action(C::ObsConfig; kwargs...)
Return the action for β = 8/g², defined as -2/g² times the sum over all traces of plaquettes.
"""
action(L::Lattice{D}, inds::Indices{D}, β::Real; kwargs...) where D = _action(plaquettesum(L, inds), β; kwargs...)
action(T::NamedTuple, β::Real; kwargs...) = action(T.lattice, T.inds, β; kwargs...)
action(C::ObsConfig; kwargs...) = _action(C.plaquette_sum, C.β; kwargs...)

function _action(plaquette_sum::Real, β::Real; log = false, iter = missing)
	log && @info "Measuring action..." iter	
	-β * plaquette_sum / 4 # 4 is due to trace normalization 
end

"""
	actionsquared(L::Lattice{D}, inds::Indices{D}, β; log = false, iter = missing) where D
	actionsquared(T::NamedTuple, β; kwargs...)
	actionsquared(C::ObsConfig; kwargs...)
Return the action squared for β = 8/g², defined as -2/g² times the sum over all traces of plaquettes.
"""
actionsquared(L::Lattice{D}, inds::Indices{D}, β::Real; kwargs...) where D = _actionsquared(plaquettesum(L, inds), β; kwargs...)
actionsquared(T::NamedTuple, β::Real; kwargs...) = actionsquared(T.lattice, T.inds, β; kwargs...)
actionsquared(C::ObsConfig; kwargs...) = _actionsquared(C.plaquette_sum, C.β; kwargs...)

function _actionsquared(plaquette_sum::Real, β::Real; log = false, iter = missing)
	log && @info "Measuring action squared..." iter
	_action(plaquette_sum, β)^2
end

#* ===== susceptibility =====
"""
	susceptibility(L::Lattice; log = false, iter = missing)
	susceptibility(T::NamedTuple; kwargs...)
	susceptibility(C::ObsConfig; kwargs...)
Return the susceptibility χ of the lattice `L`, defined as \
``\\chi = \\sum_{\\vec{x}} \\left< \\varphi_{\\vec{0}} \\right> \\left< \\varphi_{\\vec{x}} \\right> = L^D \\left< \\varphi^2 \\right>``.

Can be called with the symbol `χ`.
"""
susceptibility(L::Lattice; kwargs...) = _susceptibility(all_polyloops(L), spatialvolume(L))
susceptibility(T::NamedTuple; kwargs...) = susceptibility(T.lattice; kwargs...)
susceptibility(C::ObsConfig; kwargs...) = _susceptibility(C.polyloops, spatialvolume(C.L); kwargs...)

function _susceptibility(polyloops::Array{Float64}, volume::Int; log = false, iter = missing)
	log && @info "Measuring susceptibility..." iter
	sum(polyloops .^ 2) - sum(abs.(polyloops))^2 / volume # subtracting the condensate  
end

"""
See `susceptibility`.
"""
χ = susceptibility


"""
	susceptibility_pervolume(L::Lattice; log = false, iter = missing)
	susceptibility_pervolume(T::NamedTuple; kwargs...)
	susceptibility_pervolume(C::ObsConfig; kwargs...)
Return the susceptibility χ of the lattice `L` divided by the lattice spatial volume.

Can be called with the symbol `χᵥ`.
"""
susceptibility_pervolume(L::Lattice; kwargs...) = _susceptibility_pervolume(all_polyloops(L), spatialvolume(L))
susceptibility_pervolume(T::NamedTuple; kwargs...) = susceptibility_pervolume(T.lattice; kwargs...)
susceptibility_pervolume(C::ObsConfig; kwargs...) = _susceptibility_pervolume(C.polyloops, spatialvolume(C.L); kwargs...)

function _susceptibility_pervolume(polyloops::Array{Float64}, volume::Int; log = false, iter = missing)
	log && @info "Measuring susceptibility per volume..." iter
	_susceptibility(polyloops, volume) / volume
end

"""
See `susceptibility_pervolume`.
"""
χᵥ = susceptibility_pervolume


#* ===== Binder cumulant =====
"""
	binder(L::Lattice; log = false, iter = missing)
	binder(T::NamedTuple; kwargs...)
	binder(C::ObsConfig; kwargs...)
Return the Binder cumulant of the lattice `L`, defined as \
``g_{R} = \\frac{\\left< \\varphi^4 \\right>}{\\left< \\varphi^2 \\right>^2} - 3``.

Can be called with the symbol `gᵣ`.
"""
binder(L::Lattice; kwargs...) = _binder(all_polyloops(L), spatialvolume(L); kwargs...)
binder(T::NamedTuple; kwargs...) = binder(T.lattice; kwargs...)
binder(C::ObsConfig; kwargs...) = _binder(C.polyloops, spatialvolume(C.L); kwargs...)

function _binder(polyloops::Array{Float64}, volume::Int; log = false, iter = missing)
	log && @info "Measuring Binder cumulant..." iter
	φ⁴ = sum(polyloops .^ 4)
	φ² = sum(polyloops .^ 2)	
	volume * φ⁴ / (φ²^2) - 3
end

"""
See `binder`.
"""
gᵣ = binder


#* ===== exp. val. of Polyakov loop =====
"""
	expval_polyloop(L::Lattice; log = false, iter = missing)
	expval_polyloop(T::NamedTuple; kwargs...)
	expval_polyloop(C::ObsConfig; kwargs...)
Return the expectation value of the Polyakov loop of the lattice `L`.
"""
expval_polyloop(L::Lattice; kwargs...) = _expval_polyloop(all_polyloops(L), spatialvolume(L); kwargs...)
expval_polyloop(T::NamedTuple; kwargs...) = expval_polyloop(T.lattice; kwargs...)
expval_polyloop(C::ObsConfig; kwargs...) = _expval_polyloop(C.polyloops, spatialvolume(C.L); kwargs...)

function _expval_polyloop(polyloops::Array{Float64}, volume::Int; log = false, iter = missing) where D
	log && @info "Measuring exp. value of Polyakov loops..." iter
	sum(polyloops) / volume # sum of the traces of all polyakov loops
end

"""
	expval_polyloop_squared(L::Lattice; log = false, iter = missing)
	expval_polyloop_squared(T::NamedTuple; kwargs...)
	expval_polyloop_squared(C::ObsConfig; kwargs...)
Return the expectation value of the Polyakov loop squared of the lattice `L`.
"""
expval_polyloop_squared(L::Lattice; kwargs...) = _expval_polyloop_squared(all_polyloops(L), spatialvolume(L); kwargs...)
expval_polyloop_squared(T::NamedTuple; kwargs...) = expval_polyloop_squared(T.lattice; kwargs...)
expval_polyloop_squared(C::ObsConfig; kwargs...) = _expval_polyloop_squared(C.polyloops, spatialvolume(C.L); kwargs...)

function _expval_polyloop_squared(polyloops::Array{Float64}, volume::Int; log = false, iter = missing) where D
	log && @info "Measuring exp. value of the squares of Polyakov loops..." iter
	sum(polyloops .^ 2) / volume # sum of the traces of all polyakov loops
end

"""
	expval_modpolyloop(L::Lattice; log = false, iter = missing)
	expval_modpolyloop(T::NamedTuple; kwargs...)
	expval_modpolyloop(C::ObsConfig; kwargs...)
Return the expectation value of the modulus the Polyakov loop of the lattice `L`.
"""
expval_modpolyloop(L::Lattice; kwargs...) = _expval_modpolyloop(all_polyloops(L), spatialvolume(L); kwargs...)
expval_modpolyloop(T::NamedTuple; kwargs...) = expval_modpolyloop(T.lattice; kwargs...)
expval_modpolyloop(C::ObsConfig; kwargs...) = _expval_modpolyloop(C.polyloops, spatialvolume(C.L); kwargs...)

function _expval_modpolyloop(polyloops::Array{Float64}, volume::Int; log = false, iter = missing) where D
	log && @info "Measuring exp. value of modulus of Polyakov loops..." iter
	sum(abs.(polyloops)) / volume # sum of the traces of all polyakov loops
end


#* ===== Polyakov loops =====
"""
	corr_polyloop(L::Lattice; log = false, iter = missing)
	corr_polyloop(T::NamedTuple; kwargs...)
	corr_polyloop(C::ObsConfig; kwargs...)
Return the expectation value of the modulus the Polyakov loop of the lattice `L`.
"""
corr_polyloop(R::Int, L::Lattice; kwargs...) = _corr_polyloop(R, all_polyloops(L); kwargs...)
corr_polyloop(R::Int, T::NamedTuple; kwargs...) = corr_polyloop(R, T.lattice; kwargs...)
corr_polyloop(R::Int, C::ObsConfig; kwargs...) = _corr_polyloop(R, C.polyloops; kwargs...)

function _corr_polyloop(R::Int, polyloops::Array{Float64}; log = false, iter = missing)
	log && @info "Measuring correlation function of Polyakov loops at distance $R..." iter
	D = ndims(polyloops)
	dims = size(polyloops)
	tot = 0.0
	for x in CartesianIndices(polyloops)
		px = polyloops[x]
		for d in 1:D
			y = CartesianIndex(ntuple(i -> i == d ? mod1(x[d]+R, dims[d]) : x[d]))
			tot += px * polyloops[y]
		end
	end
	tot / (length(polyloops) * D)
end

"""
	polyloop(L::Lattice{D}, x::NTuple{Dm1, Int}; log = false, iter = missing) where {D, Dm1}
	polyloop(L, x; kwargs...)
	polyloop(L, x...; kwargs...)
	polyloop(T::NamedTuple, x; kwargs...)
	polyloop(C::ObsConfig, x...; kwargs...) 
Return the Polyakov loop calculated at spatial point `x` of the lattice `L`. `x` must be compatible with `L`'s spatial dimensions.
"""
function polyloop(L::Lattice{D}, inds::Indices{D}, x::NTuple{Dm1, Int}; log = false, iter = missing) where {D, Dm1}
	Dm1 ≠ D-1 && throw(TypeError(corr_loop, CartesianIndex{D-1}, CartesianIndex{Dm1}))

	log && @info "Measuring Polyakov loop at $x..." iter

	current_workers = vec(procs(L)) # vector of workers that own L
	dist = size(procs(L)) # how are partitions distributed among workers
	X = distribute(Array{Union{Sp2, Missing}}(undef, dist...), procs = current_workers, dist = collect(dist))
	
	with_workers(procs = current_workers) do
		if all(x .∈ tail(localindices(L)))
			time_direction = [l[1] for l in L[:L]]
			local_x = tail(Tuple(findfirst(i -> tail(Tuple(i)) == x, inds[:L])))
			X[:L][begin] = prod(time_direction[:, local_x...])
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
polyloop(L, inds, x; kwargs...) = polyloop(L, inds, Tuple(x); kwargs...)
polyloop(L, inds, x...; kwargs...) = polyloop(L, inds, Tuple(x); kwargs...)
polyloop(T::NamedTuple, x...; kwargs...) = polyloop(T.lattice, T.inds, x...; kwargs...)
polyloop(C::ObsConfig, x...; kwargs...) = polyloop(C.L, C.inds, x...; kwargs...)

"""
	twopoints_polyloop(L::Lattice{D}, x::NTuple{Dm1, Int}, y::NTuple{Dm1, Int}; log = false, iter = missing) where {D, Dm1}
	twopoints_polyloop(L, x, y; kwargs...)
	twopoints_polyloop(T::NamedTuple, x, y; kwargs...)
	twopoints_polyloop(C::ObsConfig, x, y; kwargs...)
Return the two point correlation function of the Polyakov loops at points `x` and `y`, which must be compatible with the spatial dimensions of the lattice `L`.
"""
function twopoints_polyloop(L::Lattice{D}, inds::Indices{D}, x::NTuple{Dm1, Int}, y::NTuple{Dm1, Int}; log = false, iter = missing) where {D, Dm1}
	Dm1 ≠ D-1 && throw(TypeError(corr_loop, CartesianIndex{D-1}, CartesianIndex{Dm1}))

	log && @info "Measuring two points function of Polyakov loops at $x and $y..." iter

	current_workers = vec(procs(L)) # vector of workers that own L
	dist = size(procs(L)) # how are partitions distributed among workers
	X = distribute(Array{Union{Sp2, Missing}}(undef, dist...), procs = current_workers, dist = collect(dist))
	Y = distribute(Array{Union{Sp2, Missing}}(undef, dist...), procs = current_workers, dist = collect(dist))
	
	with_workers(procs = current_workers) do
		time_direction = [l[1] for l in L[:L]]
		if all(x .∈ tail(localindices(L)))
			local_x = tail(Tuple(findfirst(i -> tail(Tuple(i)) == x, inds[:L])))
			X[:L][begin] = prod(time_direction[:, local_x...])
		end
		if all(y .∈ tail(localindices(L)))
			local_y = tail(Tuple(findfirst(i -> tail(Tuple(i)) == y, inds[:L])))
			Y[:L][begin] = prod(time_direction[:, local_y...])
		end
	end
	locX = convert(Array, X)
	locY = convert(Array, Y)
	spacedist = tail(dist)

	# only one spatial index will have a time slice of Sp2 matrices that are not missing
	trX = 0.0
	trY = 0.0
	for x in CartesianIndices(spacedist)
		if !ismissing(locX[1, x])
			trX = tr(reduce(*, locX[:, x]))
		end
		if !ismissing(locY[1, x])
			trY = tr(reduce(*, locY[:, x]))
		end
	end

	trX * trY
end
twopoints_polyloop(L, inds, x, y; kwargs...) = twopoints_polyloop(L, inds, Tuple(x), Tuple(y); kwargs...)
twopoints_polyloop(T::NamedTuple, x, y; kwargs...) = twopoints_polyloop(T.lattice, T.inds, x, y; kwargs...)
twopoints_polyloop(C::ObsConfig, x, y; kwargs...) = twopoints_polyloop(C.L, C.inds, x, y; kwargs...)

#* ===== intermediate useful functions =====
"""
	spatialvolume(L::Lattice)
Return the volume of the spatial portion of the lattice (which corresponds to all dimensions except the first one which is time).
"""
spatialvolume(L::Lattice) = prod(tail(size(L)))

"""
	plaquettesum(L::Lattice{D}, inds::Indices{D}) where D
	plaquettesum(T::NamedTuple; kwargs...)
Return the sum over the traces of all plaquettes. 
"""
function plaquettesum(L::Lattice{D}, inds::Indices{D}) where D
	current_workers = vec(procs(L))
	n_current_workers = length(current_workers)
	# initialize a vector with current_workers elements, of which each worker owns one cell
	partial = dzeros((n_current_workers,), current_workers, [n_current_workers])
	with_workers(procs = current_workers) do
		tot = 0.0
		for x in CartesianIndices(L[:L]), u in 1:D-1, v in u+1:D
			tot += tr(L[:L][x][u] * staple(L, v, u, inds[:L][x]))
		end
		partial[:L][begin] = tot
	end
	sum(partial)
end
plaquettesum(T::NamedTuple; kwargs...) = plaquettesum(T.lattice, T.inds; kwargs...)

"""
	all_polyloops(L::Lattice{D}) where D
	all_polyloops(T::NamedTuple; kwargs...)
Return an array with the same spatial dimension of the lattice `L`. Every cell contains the trace of the Polyakov loop at that spatial point.
"""
function all_polyloops(L::Lattice{D}) where D
	current_workers = vec(procs(L)) # vector of workers that own L
	dist = size(procs(L)) # how are partitions distributed among workers

	# initialize an empty array of array: each cell belongs to a worker, and it contains a D-1 dimensional array to store partial product of time-like link 
	partial = DArray(I -> Array{Array{Sp2, D-1}}(undef, length.(I)...), dist, current_workers, dist);
	
	with_workers(procs = current_workers) do
		localspacedims = tail(size(L[:L])) # size of the space part of local L
		# product along time direction for each space point
		time_direction = [l[1] for l in L[:L]] # links pointing in time direction
		partial[:L][begin] = [prod(time_direction[:, x]) for x in CartesianIndices(localspacedims)] # iterate through spatial indices of L[:L]
	end
	loc = convert(Array, partial) # bring partial into local process
	spacedist = tail(dist) # distrbutions of the space indices
	
	# product of each matrix of space points along the time axes. This is done because the time axes could be 
	# divided into different processes, and we need to multiply every partial polyakov loop.
	# list comprehension will be an array of arrays, one for each process, and each of these arrays
	# contains the product of time-link link in each lattice's position belonging to that process.
	# "mortar" transforms the array of arrays into a single array that corresponds to the spatial slice of L
	res = convert(Array, mortar([reduce(.*, loc[:, x]) for x in CartesianIndices(spacedist)]))
	tr.(res) ./ 4 # trace normalization
end
all_polyloops(T::NamedTuple; kwargs...) = all_polyloops(T.lattice; kwargs...)



#* ===== LEGACY FUNCTIONS (they are probably slower for large lattices) =====

#=

function susceptibility2(L::Lattice{D}; log = false, iter = missing) where D
	log && @info "Measuring susceptibility..." iter
	loops = all_polyloops(L)
	φ₀ = loops[begin]
	sum([φ₀ * loop for loop in loops])  
end
susceptibility2(T::NamedTuple; kwargs...) = susceptibility2(T.lattice; kwargs...)

function susceptibility_pervolume2(L::Lattice{D}; log = false, iter = missing) where D
	log && @info "Measuring susceptibility per volume..." iter
	susceptibility2(L) / spatialvolume(L)
end
susceptibility_pervolume2(T::NamedTuple; kwargs...) = susceptibility_pervolume2(T.lattice; kwargs...)
=#

#= function polyloop2(L::Lattice{D}, x::NTuple{Dm1, Int}; log = false, iter = missing) where {D, Dm1}
	Dm1 ≠ D-1 && throw(TypeError(polyloop, NTuple{D-1, Int}, NTuple{Dm1, Int}))

	log && @info "Measuring Polyakov loop at $x..." iter

	all_polyloops(L)[x...]
end
polyloop2(L, x; kwargs...) = polyloop2(L, Tuple(x); kwargs...)
polyloop2(L, x...; kwargs...) = polyloop2(L, Tuple(x); kwargs...)
polyloop2(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, x; kwargs...) = polyloop2(T[1], x; kwargs...)
polyloop2(T::NamedTuple, x; kwargs...) = polyloop2(T.lattice, x; kwargs...) =#

#= function twopoints_polyloop2(L::Lattice{D}, x::NTuple{Dm1, Int}, y::NTuple{Dm1, Int}; log = false, iter = missing) where {D, Dm1}
	Dm1 ≠ D-1 && throw(TypeError(polyloop, NTuple{D-1, Int}, NTuple{Dm1, Int}))

	log && @info "Measuring two point correlation function of Polyakov loops at $x and $y..." iter

	loops = all_polyloops(L)
	loops[x...] * loops[y...]
end
twopoints_polyloop2(L, x, y; kwargs...) = twopoints_polyloop2(L, Tuple(x), Tuple(y); kwargs...)
twopoints_polyloop2(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, x, y; kwargs...) = twopoints_polyloop2(T[1], x, y; kwargs...)
twopoints_polyloop2(T::NamedTuple, x, y; kwargs...) = twopoints_polyloop2(T.lattice, x, y; kwargs...) =#

#= function twopoints_polyloop3(L::Lattice{D}, x::NTuple{Dm1, Int}, y::NTuple{Dm1, Int}; log = false) where {D, Dm1}
	Dm1 ≠ D-1 && throw(TypeError(corr_loop, CartesianIndex{D-1}, CartesianIndex{Dm1}))

	log && @info "Measuring Polyakov loop at $x..."

	polyloop2(L, x) * polyloop2(L, y)
end
twopoints_polyloop3(L, x, y; kwargs...) = twopoints_polyloop3(L, Tuple(x), Tuple(y); kwargs...)
twopoints_polyloop3(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, x, y; kwargs...) = twopoints_polyloop3(T[1], x, y; kwargs...)
twopoints_polyloop3(T::NamedTuple, x, y; kwargs...) = twopoints_polyloop3(T.lattice, x, y; kwargs...) =#
