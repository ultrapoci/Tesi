using LinearAlgebra: I

export Link

struct Link{D} <: AbstractMatrix{ComplexF64}
	s::Sp2Element
	direction::Integer
	position::CartesianIndex{D}

	# TODO: cold and hot start
	function Link(S::Sp2Element, direction::Integer, position::CartesianIndex{D}) where D
		if abs(direction) > D 
			throw(ArgumentError("abs(direction) = $(abs(direction)) of link is bigger than dimensions = $D of the lattice"))
		end

		new{D}(S, direction, position)
	end

	function Link(S::Sp2Element, direction::Integer, position::Vararg{Integer, D}) where D
		Link(S, direction, Tuple(position))
	end
	
	function Link(S::Sp2Element, direction::Integer, position::NTuple{D, Integer}) where D
		Link(S, direction, CartesianIndex(position))
	end
	
	function Link(direction::Integer, position::CartesianIndex{D}) where D
		Link(Sp2Element(), direction, position)
	end

	function Link(direction::Integer, position::NTuple{D, Integer}) where D
		Link(Sp2Element(), direction, CartesianIndex(position))
	end

	function Link(direction::Integer, position::Vararg{Integer, D}) where D
		Link(Sp2Element(), direction, CartesianIndex(position))
	end
end
end