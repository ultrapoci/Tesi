using LinearAlgebra: I

export Link

struct Link{D} <: AbstractMatrix{ComplexF64}
	s::Sp2Element
	direction::Integer
	position::CartesianIndex{D}

	function Link(S::Sp2Element, direction::Integer, position::CartesianIndex{D}; modby::Union{Integer, NTuple{D, Integer}} = 0) where D
		if abs(direction) > D 
			throw(ArgumentError("abs(direction) = $(abs(direction)) of link is bigger than dimensions = $D of the lattice"))
		end

		if modby == 0
			new{D}(S, direction, position)
		else
			if any(d -> d ≤ 0, modby)
				throw(ArgumentError("modby = $modby in Link constructor: must be a positve integer or a tuple of positive integers."))
			end

			p = mod1.(Tuple(position), modby)
			new{D}(S, direction, CartesianIndex(p))
		end
	end

	function Link(S::Sp2Element, direction::Integer, position::Vararg{Integer, D}; modby::Union{Integer, NTuple{D, Integer}} = 0) where D
		Link(S, direction, Tuple(position); modby = modby)
	end
	
	function Link(S::Sp2Element, direction::Integer, position::NTuple{D, Integer}; modby::Union{Integer, NTuple{D, Integer}} = 0) where D
		Link(S, direction, CartesianIndex(position); modby = modby)
	end
	
	function Link(direction::Integer, position::CartesianIndex{D}; modby::Union{Integer, NTuple{D, Integer}} = 0) where D
		Link(Sp2Element(), direction, position; modby = modby)
	end

	function Link(direction::Integer, position::NTuple{D, Integer}; modby::Union{Integer, NTuple{D, Integer}} = 0) where D
		Link(Sp2Element(), direction, CartesianIndex(position); modby = modby)
	end

	function Link(direction::Integer, position::Vararg{Integer, D}; modby::Union{Integer, NTuple{D, Integer}} = 0) where D
		Link(Sp2Element(), direction, CartesianIndex(position); modby = modby)
	end
end

function Base.show(io::IO, ::MIME"text/plain", link::Link{D}) where D
    println(io, "position = $(Tuple(link.position))")
    println(io, "direction = $(link.direction)")
	show(io, "text/plain", link.s)
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

function Base.:+(l1::Link, l2::Link)
	l1.s + l2.s
end

function Base.:+(l::Link, S::Sp2Element)
	l.s + S
end

function Base.:+(S::Sp2Element, l::Link)
	S + l.s
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