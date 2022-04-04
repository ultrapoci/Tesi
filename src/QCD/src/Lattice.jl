import EnumX: @enumx

export Lattice, evensites, oddsites, directionindex, LatticeStart

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
struct Lattice{D} <: AbstractArray{Link{D}, D}
	lattice::Array{Vector{Link{D}}, D}
	
	function Lattice(::Val{D}, length::Integer; start::LatticeStart.T = LatticeStart.Cold) where D
		# TODO constraints on length and dimensions
		lattice = Array{Vector{Link}}(undef, ntuple(_ -> length, Val(D))) # `Val` is used for type stability
		
		for index in CartesianIndices(lattice)
			lattice[index] = if start ≠ LatticeStart.Empty
				links = Link{D}[]
				for direction in 1:D
					push!(links, Link(Tuple(index), direction))
				end
				links
			else
				Vector{Link}(undef, D)
			end
		end

		new{D}(lattice)
	end

	function Lattice(dimensions::Vararg{Integer, D}; start::LatticeStart.T = LatticeStart.Cold) where D
		# TODO constraints on negative dimensions
		if D == 0
			throw(ArgumentError("At least one dimension must be provided to Lattice constructor."))
		end

		lattice = Array{Vector{Link}}(undef, dimensions)

		for index in CartesianIndices(lattice)
			lattice[index] = if start ≠ LatticeStart.Empty
				links = Link{dimensions}[]
				for direction in 1:dimensions
					push!(links, Link(Tuple(index), direction))
				end
				links
			else
				Vector{Link}(undef, dimensions)
			end
		end

		new{dimensions}(lattice)
	end
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


#* ===== directionindex =====
"""
Returns a `CartesianIndex` containing all zeros except a one in the `direction` position. Basically, it returns the unit
vector in the given direction as a CartesianIndex. Use the version with `Val{dimensions}` or `dimensions` when possible: it is much quicker.
"""
function directionindex(::Val{dimensions}, direction::Integer)::CartesianIndex where dimensions
	CartesianIndex(
		ntuple(
			i -> i == direction ? 1 : 0,
			Val(dimensions)
		)
	)
end

function directionindex(dimensions::Integer, direction::Integer)::CartesianIndex
	directionindex(Val(dimensions), direction)
end

function directionindex(::Lattice{D}, direction::Integer)::CartesianIndex where D
	directionindex(Val(D), direction)
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