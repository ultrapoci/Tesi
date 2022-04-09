import EnumX: @enumx

export Lattice, evensites, oddsites, LatticeStart

#* ===== Lattice =====

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

	Lattice(Val(D), N)
Creates a D-dimensional lattice with side length N. 

	Lattice(x_0, x_1, x_2, ...)

Creates a lattice with given length for each dimension. At least one dimension must be given.
"""
struct Lattice{D} <: AbstractArray{Vector{Link{D}}, D}
	lattice::Array{Vector{Link{D}}, D}
	
	function Lattice(dimensions::NTuple{D, Integer}; start::LatticeStart.T = LatticeStart.Cold) where D
		if D == 0
			throw(ArgumentError("At least one dimension must be provided to Lattice constructor."))
		end

		if any(d -> d ≤ 0, dimensions)
			throw(ArgumentError("All dimensions provided to Lattice constructor must be strictly greater than zero."))
		end

		lattice = Array{Vector{Link}}(undef, dimensions)

		for index in CartesianIndices(lattice)
			lattice[index] = if start ≠ LatticeStart.Empty
				[Link(direction, index) for direction in 1:D]
			else
				Vector{Link}(undef, D)
			end
		end

		new{D}(lattice)
	end

	function Lattice(dimensions::Vararg{Integer, D}; start::LatticeStart.T = LatticeStart.Cold) where D
		Lattice(Tuple(dimensions); start = start)
	end

	function Lattice(::Val{D}, length::Integer; start::LatticeStart.T = LatticeStart.Cold) where D
		Lattice(ntuple(_ -> length, Val(D)); start = start)
	end
end

function Base.getindex(L::Lattice, i::CartesianIndex)
	t = Tuple(i)
	index = mod1.(t, size(L)) |> CartesianIndex
	getindex(L.lattice, index)
end

function Base.setindex!(L::Lattice, v, i::CartesianIndex)
	t = Tuple(i)
	index = mod1.(t, size(L)) |> CartesianIndex
	setindex!(L.lattice, v, index)
end

# Linear indexing, without `mod`
function Base.getindex(L::Lattice, i)
	getindex(L.lattice, i)
end

# Multiple indices: every entry is `mod`ed by the dimensions of the lattice to wrap around
function Base.getindex(L::Lattice, i...)
	getindex(L.lattice, mod1.(i, size(L))...)
end

# Linear indexing, without `mod`
function Base.setindex!(L::Lattice, v, i)
	setindex!(L.lattice, v, i)
end

# Multiple indices: every entry is `mod`ed by the dimensions of the lattice to wrap around
function Base.setindex!(L::Lattice, v, i...)
	setindex!(L.lattice, v, mod1.(i, size(L))...)
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

#* ===== Iteration over even and odd sites =====
"""
Iterator over a lattice's even or odd sites, where even (odd) site means that the sum over all coordinates of the site is even (odd)  
"""
struct EvenOddLattice
	lattice::Lattice
	mod2_result::Int # it is 1 if odd sites, or 0 if even sites
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
	iterate(L, Tuple(CartesianIndices(L.lattice)))
end

function Base.iterate(L::EvenOddLattice, state)
	lattice = L.lattice
	mod2_result = L.mod2_result

	indices_tail = state
	for cartesian_index in state
		indices_tail = Base.tail(indices_tail)
		s = sum(Tuple(cartesian_index)) # sum over all coordinates of site's position
		if s % 2 == mod2_result # = 0 if even, = 1 if odd
			return (lattice[cartesian_index], indices_tail)
		end
	end

	nothing	
end

function directionindex(::Lattice{D}, direction::Integer)::CartesianIndex where D
	directionindex(Val(D), direction)
end