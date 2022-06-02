const Mask{D} = DArray{Bool, D, Array{Bool, D}}
const Indices{D} = DArray{CartesianIndex{D}, D, Array{CartesianIndex{D}, D}}
const Site{D} = MVector{D, Sp2}
const LocalLattice{D} = Array{Site{D}, D}
const Lattice{D} = DArray{Site{D}, D, LocalLattice{D}}

"""
	newlattice(dims::Dims{D}; start = :cold, kwargs...) where D
	newlattice(dims::Vararg{Int, D}; kwargs...) where D
Create a D-dimensional lattice of dimensions `dims`, where each cell contains a `MVector` (mutable, static vector) of \
length `D`. Each entry of the vector contains an `Sp2` object and represents a link of the lattice pointing in a given \
direction. The first entry is the time direction, while the entries from 2 to `D` are the spatial directions. 

Returns a `NamedTuple` containing:
- `lattice`, which is a `DArray` built as explained above
- `mask`, which is a `DArray` of booleans, where `true` indicates the corresponding site is even, and `false` indicates the \
corresponding site is odd
- `inds`, which is a `DArray` of `CartesianIndex` objects, which describes the index position of the cell they belong to

`mask` and `inds` are distributed among processes as `lattice` and have the same dimensions.

Any keyword argument passed to `newlattice` is passed to the function `distribute` from DistributedArrays.jl.
"""
function newlattice(dims::Dims{D}; start = :cold, kwargs...) where D
	lattice = if start == :cold
		distribute([MVector{D}([Sp2() for _ in 1:D]) for _ in CartesianIndices(dims)]; kwargs...)
	elseif start == :hot
		distribute([MVector{D}([rand(Sp2) for _ in 1:D]) for _ in CartesianIndices(dims)]; kwargs...)
	else
		throw(ArgumentError("start must be :hot or :cold, got :$start."))
	end

	mask = distribute(evenmask(lattice), lattice)
	inds = distribute(CartesianIndices(lattice), lattice)

	(lattice = lattice, mask = mask, inds = inds)
end
newlattice(dims::Vararg{Int, D}; kwargs...) where D = newlattice(Tuple(dims); kwargs...)

evenmask(A) = [iseven(sum(Tuple(i))) for i in CartesianIndices(A)]
oddmask(A) = [isodd(sum(Tuple(i))) for i in CartesianIndices(A)]

addtuple(n::Int, d::Int, p::Dims{D}) where D = ntuple(i -> i == d ? p[i]+n : p[i], D)::Dims{D}
addtuple(n::Int, d::Int, p::Vararg{Int, D}) where D = addtuple(n, d, Tuple(p))::Dims{D}

"""
	getlink(L::Lattice{D}, u::Int, x::Dims{D}) where D
	getlink(L::Lattice{D}, u::Int, x::Vararg{Int, D}) where D
	getlink(L::Lattice{D}, u::Int, x::CartesianIndex{D}) where D
Return the link (an `Sp2` object) in the `x` site of the lattice `L` pointing in the `u` direction. 

`x` is modded by `L`'s dimensions as to wrap around if it doesn't fit the correct dimensions, simulating periodic \
boundary conditions.

`u` can be negative: in this case, the link pointing in a negative direction is the link coming from a near neighbour of \
`x`, pointing towards `x`, which is then Hermitian conjugated.

"""
function getlink(L::Lattice{D}, u::Int, x::Dims{D})::Sp2 where D
	if u ∈ 1:D
		p = CartesianIndex(mod1.(x, size(L)))
		@inbounds L[p][u]
	elseif u ∈ -D:-1
		p = CartesianIndex(mod1.(addtuple(-1, -u, x), size(L)))
		@inbounds L[p][-u]'
	else
		throw(ArgumentError("Link direction for lattice must be in range [1, $D] or [-$D, -1]. Got u = $u."))
	end
end
getlink(L::Lattice{D}, u::Int, x::Vararg{Int, D}) where D = getlink(L, u, Tuple(x))
getlink(L::Lattice{D}, u::Int, x::CartesianIndex{D}) where D = getlink(L, u, Tuple(x))

"""
	staple(L::Lattice{D}, d::Int, u::Int, x::Dims{D}) where D
Calculate the staple in the direction `v` around the link at position `x` pointing in the direction `u`.

Note that `v` ≠ `u`. `u` must be positive and `u` ≤ `D`, while `v` must be in range [1, `D`] or [-`D`, -1].
"""
function staple(L::Lattice{D}, v::Int, u::Int, x::Dims{D}) where D
	if u ∉ 1:D
		throw(ArgumentError("The link's direction must be in range [1, $D]: got u=$u."))
	elseif v ∉ 1:D && v ∉ -D:-1
		throw(ArgumentError("The staple's direction must be in range [1, $D] or [-$D, -1], got d=$v."))
	elseif u == v
		throw(ArgumentError("The link's direction (u=$v) cannot be the same as the staple's direction."))
	end

	û = MVector(x)
	v̂ = copy(û)
	ŵ = copy(û) # û + v̂
	
	sv = sign(v)
	iv = abs(v)
	
	@inbounds û[u] += 1
	@inbounds v̂[iv] += sv
	@inbounds ŵ[u] += 1
	@inbounds ŵ[iv] += sv

	U₁ = getlink(L,  v, Tuple(û))
	U₂ = getlink(L, -u, Tuple(ŵ))
	U₃ = getlink(L, -v, Tuple(v̂))

	U₁ * U₂ * U₃
end
staple(L::Lattice{D}, v::Int, u::Int, x::Vararg{Int, D}) where D = staple(L, v, u, Tuple(x))
staple(L::Lattice{D}, v::Int, u::Int, x::CartesianIndex{D}) where D = staple(L, v, u, Tuple(x))

"""
	sumstaples(L::Lattice{D}, u::Int, x::Dims{D}) where D
	sumstaples(L::Lattice{D}, u::Int, x::Vararg{Int, D}) where D
	sumstaples(L::Lattice{D}, u::Int, x::CartesianIndex{D}) where D
Returns the sum of all the staples surrounding the link at position `x` of the lattice `L`, pointing in the `u` direction.
`u` must be positive. 
"""
function sumstaples(L::Lattice{D}, u::Int, x::Dims{D})::SMatrix{4, 4, ComplexF64} where D
	u ∉ 1:D && throw(ArgumentError("Link's direction u must be in range [1, D], got $u."))	

	r = mod1.(u+1:u+D-1, D) # generate all D directions except 'u'
	sum(staple(L, v, u, x) for v in vcat(r, -r))
end
sumstaples(L::Lattice{D}, u::Int, x::Vararg{Int, D}) where D = sumstaples(L, u, Tuple(x))
sumstaples(L::Lattice{D}, u::Int, x::CartesianIndex{D}) where D = sumstaples(L, u, Tuple(x))

plaquette(L::Lattice{D}, v::Int, u::Int, x::Dims{D}) where D = getlink(L, u, x) * staple(L, v, u, x)
plaquette(L::Lattice{D}, v::Int, u::Int, x::Vararg{Int, D}) where D = plaquette(L, v, u, Tuple(x))
plaquette(L::Lattice{D}, v::Int, u::Int, x::CartesianIndex{D}) where D = plaquette(L, v, u, Tuple(x))

#=
function sumstaples(L::Lattice{D}, u::Int, x::Dims{D})::SMatrix{4, 4, ComplexF64} where D
	u ∉ 1:D && throw(ArgumentError("Link's direction u must be in range [1, D], got $u."))	
	
	total = zeros(MMatrix{4, 4, ComplexF64})
	for v in mod1.(u+1:u+D-1, D) # generate all D directions except u
		total += staple(L, v, u, x) + staple(L, -v, u, x)
	end

	SMatrix(total)
end
sumstaples(L::Lattice{D}, u::Int, x::Vararg{Int, D}) where D = sumstaples(L, u, Tuple(x))
sumstaples(L::Lattice{D}, u::Int, x::CartesianIndex{D}) where D = sumstaples(L, u, Tuple(x))
=#