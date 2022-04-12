export Lattice, evensites, oddsites, evenlinks, oddlinks, updatelattice!

#* ===== Lattice =====
"""
Creates a D-dimensional lattice with side length N.

	Lattice(Val(D), N[; start::Symbol = :cold])

Creates a lattice with given length for each dimension. At least one dimension must be given.

	Lattice(x_0, x_1, x_2, ... [; start::Symbol = :cold])
	Lattice((x_0, x_1, x_2)[; start::Symbol = :cold])

If `start` = :cold, all links of the lattice are set to be the identity element.
If `start` = :hot, all links are randomized and normalized Sp2Element.
If `start` = :empty, the lattice contains initialized vectors of undef links.

"""
struct Lattice{D} <: AbstractArray{Vector{Link{D}}, D}
	lattice::Array{Vector{Link{D}}, D}
	
	function Lattice(dimensions::NTuple{D, Integer}; start::Symbol = :cold) where D
		if start ∉ [:cold, :hot, :empty]
			throw(ArgumentError("start keyword argument must be equal to :empty, :cold or :hot. Got start = $start."))
		elseif D == 0
			throw(ArgumentError("At least one dimension must be provided to Lattice constructor."))
		elseif any(d -> d ≤ 0, dimensions)
			throw(ArgumentError("All dimensions provided to Lattice constructor must be strictly greater than zero."))
		end

		lattice = Array{Vector{Link}}(undef, dimensions)

		for index in CartesianIndices(lattice)
			lattice[index] = if start == :cold
				[Link(direction, index) for direction in 1:D]
			elseif start == :hot
				[Link(randSp2(), direction, index) for direction in 1:D]
			else
				Vector{Link}(undef, D)
			end			
		end

		new{D}(lattice)
	end

	function Lattice(dimensions::Vararg{Integer, D}; start::Symbol = :cold) where D
		Lattice(Tuple(dimensions); start = start)
	end

	function Lattice(::Val{D}, length::Integer; start::Symbol = :cold) where D
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

#* ===== updatelattice! =====
function updatelattice!(L::Lattice{D}, link::Link{D}) where D
	L[link.position][link.direction] = link
end

function updatelattice!(L::Lattice{D}, links::Vector{Link{D}}) where D
	for link in links
		updatelattice!(L, link)
	end
end

#* ===== Iteration over even and odd sites and links =====
"""
Iterator over a lattice's even or odd sites, where even (odd) site means that the sum over all coordinates of the site is even (odd).
"""
struct EvenOddLattice{D}
	lattice::Lattice{D}
	mod2_result::Int # it is 1 if odd sites, or 0 if even sites

	function EvenOddLattice(lattice::Lattice{D}, mod2_result::Int) where D
		if mod2_result ≠ 0 && mod2_result ≠ 1
			throw(ArgumentError("mod2_result in EvenOddLattice must be 0 or 1, got $mod2_result."))
		end
		new{D}(lattice, mod2_result)
	end
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

struct EvenOddLinks{D}
	eolattice::EvenOddLattice{D}
	direction::Integer

	function EvenOddLinks(eolattice::EvenOddLattice{D}, direction::Integer) where D
		if direction ∉ 1:D
			throw(ArgumentError("Direction given is $direction, but it must be in range [1, $D]."))
		end
		new{D}(eolattice, direction)
	end
end

evenlinks(L::Lattice, direction::Integer) = EvenOddLinks(evensites(L), direction)
oddlinks(L::Lattice, direction::Integer) = EvenOddLinks(oddsites(L), direction)

function Base.iterate(L::EvenOddLinks)
	iterate(L, Tuple(CartesianIndices(L.eolattice.lattice)))
end

function Base.iterate(L::EvenOddLinks, state)
	result = iterate(L.eolattice, state)
	if !isnothing(result)
		return (result[1][L.direction], result[2])
	end
	nothing
end