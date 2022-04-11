module QCD

include("SU2Element.jl")
include("Sp2Element.jl")
include("Link.jl")
include("Lattice.jl")

export getlink, staple, plaquette, unittuple, directionindex

function unittuple(::Val{dimensions}, direction::Integer) where dimensions
	if dimensions ≤ 0
		throw(ArgumentError("dimensions = $dimensions in unittuple must be a positive integer"))
	end
	
	ntuple(
		i -> i == abs(direction) ? sign(direction) : 0,
		Val(dimensions)
	)
end

function unittuple(dimensions::Integer, direction::Integer)
	unittuple(Val{dimensions}, direction)
end


#* ===== directionindex =====
"""
	directionindex(::Val{dimensions}, direction::Integer)::CartesianIndex where dimensions
	directionindex(dimensions::Integer, direction::Integer)::CartesianIndex
	directionindex(::Lattice{D}, direction::Integer)::CartesianIndex where D

Returns a `CartesianIndex` containing all zeros except a one in the `direction` position. Basically, it returns the unit
vector in the given direction as a CartesianIndex. Use the version with `Val{dimensions}` or `dimensions` when possible: it is much quicker.
It accounts for negative directions: if `direction` is < 0, it will return an index with a -1 instead of a 1 in `abs(direction)` position.
"""
function directionindex(::Val{dimensions}, direction::Integer)::CartesianIndex where dimensions
	CartesianIndex(unittuple(Val(dimensions), direction))
end

function directionindex(dimensions::Integer, direction::Integer)::CartesianIndex
	directionindex(Val(dimensions), direction)
end

function directionindex(::Lattice{D}, direction::Integer)::CartesianIndex where D
	directionindex(Val(D), direction)
end


#* ===== getlink =====
"""
	getlink(lattice::Lattice{D}, direction::Integer, position::CartesianIndex{D}) where D
Returns the link at the given `position` in the given `direction`. If the `direction` is negative, 
returns the adjoint of the link that points towards the given `position` from the neighbour. 
"""
function getlink(lattice::Lattice{D}, direction::Integer, position::CartesianIndex{D}) where D
	if direction ∈ 1:D 
		lattice[position][direction]
	elseif direction ∈ -D:-1
		p = position + directionindex(D, direction)
		s = lattice[p][-direction].s
		Link(s', direction, position; modby = size(lattice))
	else
		throw(ArgumentError("direction = $direction must be in the interval [1, $D] or [-$D, -1]."))
	end
end

function getlink(lattice::Lattice{D}, direction::Integer, position::NTuple{D, Integer}) where D
	getlink(lattice, direction, CartesianIndex(position))
end

function getlink(lattice::Lattice{D}, direction::Integer, position::Vararg{Integer, D}) where D
	getlink(lattice, direction, CartesianIndex(position))
end


#* ===== staple =====
"""
	staple(lattice::Lattice{D}, link::Link{D}, direction::Integer) where D

Returns the staple around `link` in the given `direction`. If, for example, the link points up (↑), 
and `direction` is positive, this returns:

→ * ↓ * ←

The returned product already follows the given link direction, meaning that multiplying this by the link itself
gives the complete plaquette. Only accepts Link with positive direction, and `direction` must not be equal to the
Link's direction.
"""
function staple(lattice::Lattice{D}, link::Link{D}, direction::Integer) where D
	if link.direction ≤ 0 
		throw(ArgumentError("The provided link has negative direction = $(link.direction)."))
	elseif direction ∉ 1:D && direction ∉ -D:-1
		throw(ArgumentError("The staple's direction = $direction must be in range [1, $D] or [-$D, -1]."))
	elseif link.direction == direction
		throw(ArgumentError("Given direction = $direction cannot be the same as the link's direction."))
	end

	x = link.position
	u = link.direction
	v = direction

	û = directionindex(D, u)
	v̂ = directionindex(D, v)
	
	U₁ = getlink(lattice,  v, x + û).s
	U₂ = getlink(lattice, -u, x + û + v̂).s
	U₃ = getlink(lattice, -v, x + v̂).s

	U₁ * U₂ * U₃
end

#* ===== plaquette =====

function plaquette(lattice::Lattice{D}, link::Link{D}, direction::Integer) where D
	#= if link.direction ≤ 0 
		throw(ArgumentError("The provided link has negative direction = $(link.direction)."))
	elseif direction ∉ 1:D && direction ∉ -D:-1
		throw(ArgumentError("The plaquette's direction = $direction must be in range [1, $D] or [-$D, -1]."))
	elseif link.direction == direction
		throw(ArgumentError("Given direction = $direction cannot be the same as the link's direction."))
	end =#

	R = staple(lattice, link, direction)
	link.s * R
end

end # module QCD