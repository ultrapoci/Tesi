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

function Base.:*(l1::Link, l2::Link)
	l1.s * l2.s
end

function Base.:*(l::Link, S::Sp2Element)
	l.s * S
end

function Base.:*(S::Sp2Element, l::Link)
	S * l.s
end

function Base.getindex(link::Link, i...)
	getindex(link.s, i...)
end

function Base.setindex!(link::Link, v, i...)
	setindex!(link.s, v, i...)
end

function Base.firstindex(link::Link)
	firstindex(link.s)
end

function Base.lastindex(link::Link)
	lastindex(link.s)
end

function Base.size(link::Link)
	size(link.s)
end