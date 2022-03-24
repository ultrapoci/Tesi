module Tesi

using DrWatson
DrWatson.@quickactivate "Tesi" 

import LinearAlgebra: I
import EnumX: @enumx

export Link, Lattice, evensites, oddsites, directionindex, LatticeStart


# ==== Link =====

struct Link{D}
	position::NTuple{D, Integer}
	direction::Integer
	m::Matrix{Complex} # TODO: see if there's a better way to store Sp(2) in memory

	# TODO: cold and hot start
	function Link(position::NTuple{D, Integer}, direction::Integer) where D
		if direction > D 
			throw(ArgumentError("direction = $direction of given link is bigger than dimensions = $D of the lattice"))
		end
		m = Matrix{Complex}(I, 4, 4) # TODO: see if Hermitian(m) is useful 
		new{D}(position, direction, m)
	end            
end


# =====  Lattice =====

@enumx LatticeStart begin 
	"""
	Creates a lattice as an `Array{Vector{Link{D}}, D}`, where every `Vector{Link{D}}` is initialized with `undef`.
	"""
	Empty

	"""
	Creates a lattice where all `Link`s are initialized with the identity matrix.
	"""
	Cold

	"""
	Creates a lattice where all `Link`s are initialized with a random matrix.
	"""
	Hot
end

"""
Data structure containing the lattice as `dimensions`-dimensional array.
"""
struct Lattice{D} <: AbstractArray{Link{D}, D}
	dimensions::Integer
	length::Integer
	lattice::Array{Vector{Link{D}}, D}
	
	#function Lattice(dimensions::Integer, length::Integer, start::LatticeStart.T = LatticeStart.Cold) 
	function Lattice(::Val{dimensions}, length::Integer, start::LatticeStart.T = LatticeStart.Cold) where dimensions
		# TODO constraints on length and dimensions
		lattice = Array{Vector{Link}}(undef, ntuple(_ -> length, Val(dimensions))...) # `Val` is used for type stability
		
		for index in CartesianIndices(lattice)
			lattice[index] = if start ≠ LatticeStart.Empty
				links = Link[]
				for direction in 1:dimensions
					push!(links, Link(Tuple(index), direction))
				end
				links
			else
				Vector{Link}(undef, dimensions)
			end
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

"""
Returns a `CartesianIndex` containing all zeros except a one in the `direction` position. Basically, it returns the unit
vector in the given direction as a CartesianIndex. Use the version with `Val{dimensions}` when possible: it is much quicker.
"""
function directionindex(L::Lattice, direction::Integer)
	directionindex(L.dimensions, direction)
end

#function directionindex(dimensions::Integer, direction::Integer)
function directionindex(::Val{dimensions}, direction::Integer) where dimensions
	CartesianIndex(
		ntuple(
			i -> i == direction ? 1 : 0,
			Val(dimensions)
		)
	)
end

# ===== Iteration over even and odd sites =====
"""
Iterator over a lattice's even or odd sites, where even (odd) site means that the sum over all coordinates of the site is even (odd)  
"""
struct EvenOddLattice
	lattice::Lattice
	mod_result::Int # it is 1 if odd sites, or 0 if even sites
end

"""
Returns an iterator over all even sites of a given lattice. A site is even if the sum of all its coordinates are even. 
This is useful in conjuction with `oddsites`: even and odd sites don't influence one another when applying the updating algorithm. 
"""
evensites(L::Lattice) = EvenOddLattice(L, 0)

"""
Returns an iterator over all odd sites of a given lattice. A site is odd if the sum of all its coordinates are odd. 
This is useful in conjuction with `evensites`: even and odd sites don't influence one another when applying the updating algorithm.
"""
oddsites(L::Lattice) = EvenOddLattice(L, 1)

function Base.iterate(L::EvenOddLattice)
	# decide whether the first site is even or odd
	if L.lattice.dimensions % 2 == L.mod_result # = 0 if even, = 1 if odd
		(L.lattice[begin], Base.tail(Tuple(CartesianIndices(L.lattice))))
	else
		(L.lattice[begin + 1], Base.tail(Base.tail(Tuple(CartesianIndices(L.lattice)))))
	end
end

function Base.iterate(L::EvenOddLattice, state)
	indices_tail = state
	for cartesian_index in state
		indices_tail = Base.tail(indices_tail)
		s = sum(Tuple(cartesian_index)) # sum over all coordinates of site's position
		if s % 2 == L.mod_result # = 0 if even, = 1 if odd
			return (L.lattice[cartesian_index], indices_tail)
		end
	end

	nothing	
end

end # module Tesi