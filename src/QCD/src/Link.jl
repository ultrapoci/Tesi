using LinearAlgebra: I

export Link

struct Link{D} <: AbstractMatrix{ComplexF64}
	s::Sp2Element
	direction::Integer
	position::CartesianIndex{D}

	function Link(S::Sp2Element, direction::Integer, position::CartesianIndex{D}; modby::Union{Integer, NTuple{D, Integer}} = 0) where D
		if direction ∉ 1:D
			throw(ArgumentError("direction = $direction of link is must be in range [1, $D]."))
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
	
	function Link(direction::Integer, position::CartesianIndex{D}; type::Type{<:Sp2Element} = Sp2ElementA, modby::Union{Integer, NTuple{D, Integer}} = 0) where D
		Link(Sp2Element(type), direction, position; modby = modby)
	end

	function Link(direction::Integer, position::NTuple{D, Integer}; type::Type{<:Sp2Element} = Sp2ElementA, modby::Union{Integer, NTuple{D, Integer}} = 0) where D
		Link(Sp2Element(type), direction, CartesianIndex(position); modby = modby)
	end

	function Link(direction::Integer, position::Vararg{Integer, D}; type::Type{<:Sp2Element} = Sp2ElementA, modby::Union{Integer, NTuple{D, Integer}} = 0) where D
		Link(Sp2Element(type), direction, CartesianIndex(position); modby = modby)
	end
end

function Base.show(io::IO, ::MIME"text/plain", link::Link{D}) where D
    println(io, "position = $(Tuple(link.position))")
    println(io, "direction = $(link.direction)")
	show(io, "text/plain", link.s)
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