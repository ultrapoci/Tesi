module Tesi

using DrWatson
DrWatson.@quickactivate "Tesi" 

import LinearAlgebra

export Link, Lattice, evensites, oddsites

struct Link{D}
	position::NTuple{D, Int}
	direction::Int
	m::Matrix{Complex} # TODO: see if there's a better way to store Sp(2) in memory

	function Link(position::NTuple{D, Int}, direction::Int) where D
		if direction > D 
			throw(ArgumentError("direction = $direction of given link is bigger than dimensions = $D of the lattice"))
		end
		m = Matrix{Complex}(LinearAlgebra.I, 4, 4) # TODO: see if Hermitian(m) is useful 
		new{D}(position, direction, m)
	end            
end

"""
Data structure containing the lattice as `dimensions`-dimensional array.
"""
struct Lattice{D} <: AbstractArray{Link{D}, D}
	dimensions::Int
	length::Int
	lattice::Array{Vector{Link{D}}, D}
	
	function Lattice(dimensions::Int, length::Int) # TODO: cold and hot start
		# TODO constraints on length and dimensions
		lattice = Array{Vector{Link}}(undef, tuple(fill(length, dimensions)...)) 
		for index in CartesianIndices(lattice)
			links = Link[]
			for direction in 1:dimensions
				push!(links, Link(Tuple(index), direction))
			end
			lattice[index] = links
		end
		new{dimensions}(dimensions, length, lattice)
	end
end

# Linear indexing, without `mod`
function Base.getindex(L::Lattice, i)
	getindex(L.lattice, i)
end

# Multiple indices: every entry is `mod`ed by the length of the lattice to wrap around
function Base.getindex(L::Lattice, i...)
	getindex(L.lattice, mod1.(i, L.length)...)
end

# Linear indexing, without `mod`
function Base.setindex!(L::Lattice, v, i)
	setindex!(L.lattice, v, i)
end

# Multiple indices: every entry is `mod`ed by the length of the lattice to wrap around
function Base.setindex!(L::Lattice, v, i...)
	setindex!(L.lattice, v, mod1.(i, L.length)...)
end

function Base.firstindex(L::Lattice)
	firstindex(L.lattice)
end

function Base.lastindex(L::Lattice)
	lastindex(L.lattice)
end

function Base.size(L::Lattice)
	size(L.lattice)
end


# ===== Iteration over even sites =====

struct EvenLattice
	lattice::Lattice
end

evensites(L) = EvenLattice(L)

function Base.iterate(L::EvenLattice)
	if L.lattice.dimensions % 2 == 0
		(L.lattice[begin], Base.tail(Tuple(CartesianIndices(L.lattice))))
	else
		(L.lattice[begin + 1], Base.tail(Base.tail(Tuple(CartesianIndices(L.lattice)))))
	end
end

function Base.iterate(L::EvenLattice, state)
	if length(state) == 0
		return nothing
	end

	indices_tail = state
	for cartesian_index in state
		indices_tail = Base.tail(indices_tail)
		s = sum(Tuple(cartesian_index)) # sum over all coordinates of link's position
		if s % 2 == 0 # it's even
			return (L.lattice[cartesian_index], indices_tail)
		end
	end

	nothing	
end


# ===== Iteration over odd sites =====

struct OddLattice
	lattice::Lattice
end

oddsites(L) = OddLattice(L)

function Base.iterate(L::OddLattice)
	if L.lattice.dimensions % 2 == 1 
		(L.lattice[begin], Base.tail(Tuple(CartesianIndices(L.lattice))))
	else
		(L.lattice[begin + 1], Base.tail(Base.tail(Tuple(CartesianIndices(L.lattice)))))
	end
end

function Base.iterate(L::OddLattice, state)
	if length(state) == 0
		return nothing
	end

	indices_tail = state
	for cartesian_index in state
		indices_tail = Base.tail(indices_tail)
		s = sum(Tuple(cartesian_index)) # sum over all coordinates of link's position
		if s % 2 == 1 # it's odd
			return (L.lattice[cartesian_index], indices_tail)
		end
	end

	nothing	
end

end # module Tesi