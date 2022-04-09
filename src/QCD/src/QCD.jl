module QCD

export unittuple, directionindex, getlink

function unittuple(::Val{dimensions}, direction::Integer) where dimensions
	ntuple(
		i -> i == direction ? 1 : 0,
		Val(dimensions)
	)
end

#* ===== directionindex =====
"""
Returns a `CartesianIndex` containing all zeros except a one in the `direction` position. Basically, it returns the unit
vector in the given direction as a CartesianIndex. Use the version with `Val{dimensions}` or `dimensions` when possible: it is much quicker.
"""
function directionindex(::Val{dimensions}, direction::Integer)::CartesianIndex where dimensions
	CartesianIndex(
		unittuple(Val(dimensions), direction)
	)
end

function directionindex(dimensions::Integer, direction::Integer)::CartesianIndex
	directionindex(Val(dimensions), direction)
end

include("Sp2Element.jl")
include("Link.jl")
include("Lattice.jl")

"""
Returns the link at the given position in the given direction. If the direction is negative, 
returns the adjoint of the link that points towards the given position. 
"""
function getlink(lattice::Lattice{D}, direction::Integer, position::CartesianIndex{D}) where D
	if direction ∈ 1:D 
		lattice[position][direction]
	elseif direction ∈ -D:-1
		d = -direction
		p = position - directionindex(Val(D), d)
		s = lattice[p][d].s
		Link(s', direction, position)
	else
		throw(ArgumentError("Direction = $direction must be in the interval [1, $D] or [-$D, -1]."))
	end
end

function getlink(lattice::Lattice{D}, direction::Integer, position::NTuple{D, Integer}) where D
	getlink(lattice, direction, CartesianIndex(position))
end

function getlink(lattice::Lattice{D}, direction::Integer, position::Vararg{Integer, D}) where D
	getlink(lattice, direction, CartesianIndex(position))
end

end # module QCD