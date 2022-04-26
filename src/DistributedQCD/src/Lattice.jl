const Mask{D} = DArray{Bool, D, Array{Bool, D}}
const Indices{D} = DArray{CartesianIndex{D}, D, Array{CartesianIndex{D}, D}}
const Site{D} = MVector{D, Sp2}
const LocalLattice{D} = Array{Site{D}, D}
const Lattice{D} = DArray{Site{D}, D, LocalLattice{D}}

function newlattice(dims::NTuple{D, Int}; start = :cold, kwargs...) where D
	lattice = if start == :cold
		distribute([MVector{D}([Sp2() for _ in 1:D]) for _ in CartesianIndices(dims)]; kwargs...)
	elseif start == :hot
		distribute([MVector{D}([rand(Sp2) for _ in 1:D]) for _ in CartesianIndices(dims)]; kwargs...)
	else
		throw(ArgumentError("start must be :hot or :cold, got :$start."))
	end

	mask = distribute(evenmask(lattice), lattice)
	inds = distribute(indices(lattice), lattice)

	(lattice = lattice, mask = mask, inds = inds)
end
newlattice(dims::Vararg{Int, D}; kwargs...) where D = newlattice(Tuple(dims); kwargs...)

addtuple(n::Int, d::Int, p::Tuple{Vararg{Int}}) = tuple(p[begin:d-1]..., p[d]+n, p[d+1:end]...)
addtuple(n::Int, d::Int, p::Vararg{Int}) = addtuple(n, d, Tuple(p))

function getlink(L::Lattice{D}, u::Int, x::NTuple{D, Int}) where D
	if u ∈ 1:D
		p = CartesianIndex(mod1.(x, size(L)))
		@inbounds L[p][u]
	elseif u ∈ -D:-1
		p = CartesianIndex(mod1.(addtuple(-1, -u, x), size(L)))
		@inbounds L[p][-u]'
	else
		throw(ArgumentError("Link direction for lattice must be in range [1, $D] or [-$D, -1]. Got u = $u"))
	end
end
getlink(L::Lattice{D}, u::Int, x::Vararg{Int, D}) where D = getlink(L, u, Tuple(x))
getlink(L::Lattice{D}, u::Int, x::CartesianIndex{D}) where D = getlink(L, u, Tuple(x))

evenmask(A) = [iseven(sum(Tuple(i))) for i in CartesianIndices(A)]
oddmask(A) = [isodd(sum(Tuple(i))) for i in CartesianIndices(A)]
indices(A) = [i for i in CartesianIndices(A)]

"""
	staple(L::Lattice{D}, d::Int, u::Int, x::NTuple{D, Int}) where D
Calculate the staple in the direction `d` around the link at position `x` pointing in the direction `u`.

Note that `d` ≠ `u`. `u` must be positive and `u` ≤ `D`, while `d` must be in range [1, `D`] or [-`D`, -1].
"""
function staple(L::Lattice{D}, v::Int, u::Int, x::NTuple{D, Int}) where D
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

	su = sign(u)
	sv = sign(v)
	iu = abs(u)
	iv = abs(v)

	@inbounds û[iu] += su
	@inbounds v̂[iv] += sv
	@inbounds ŵ[iu] += su
	@inbounds ŵ[iv] += sv

	U₁ = getlink(L,  v, Tuple(û))
	U₂ = getlink(L, -u, Tuple(ŵ))
	U₃ = getlink(L, -v, Tuple(v̂))

	U₁ * U₂ * U₃
end
staple(L::Lattice{D}, v::Int, u::Int, x::Vararg{Int, D}) where D = staple(L, v, u, Tuple(x))
staple(L::Lattice{D}, v::Int, u::Int, x::CartesianIndex{D}) where D = staple(L, v, u, Tuple(x))

function sumstaples(L::Lattice{D}, u::Int, x::NTuple{D, Int}) where D
	if u ∉ 1:D
		throw(ArgumentError("Link's direction u must be in range [1, D], got $u."))	
	end

	total = zeros(MMatrix{4, 4, ComplexF64})
	for v in mod1.(u+1:u+D-1, D) # generate all D directions except u
		total += staple(L, v, u, x) + staple(L, -v, u, x)
	end
	total
end
sumstaples(L::Lattice{D}, u::Int, x::Vararg{Int, D}) where D = sumstaples(L, u, Tuple(x))
sumstaples(L::Lattice{D}, u::Int, x::CartesianIndex{D}) where D = sumstaples(L, u, Tuple(x))

plaquette(L::Lattice{D}, v::Int, u::Int, x::NTuple{D, Int}) where D = getlink(L, u, x) * staple(L, v, u, x)
plaquette(L::Lattice{D}, v::Int, u::Int, x::Vararg{Int, D}) where D = plaquette(L, v, u, Tuple(x))
plaquette(L::Lattice{D}, v::Int, u::Int, x::CartesianIndex{D}) where D = plaquette(L, v, u, Tuple(x))
