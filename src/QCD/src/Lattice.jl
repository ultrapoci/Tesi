export Lattice, evensites, oddsites, evenlinks, oddlinks, iterlinks, updatelattice!, everylink

#* ===== Lattice =====
"""
Creates a D-dimensional lattice with side length N.

	Lattice(Val(D), N[; type::Type{<:Sp2Element} = Sp2ElementA, start::Symbol = :cold])

Creates a lattice with given length for each dimension. At least one dimension must be given.

	Lattice(x_0, x_1, x_2, ... [; type::Type{<:Sp2Element} = Sp2ElementA, start::Symbol = :cold])
	Lattice((x_0, x_1, x_2)[; type::Type{<:Sp2Element} = Sp2ElementA, start::Symbol = :cold])

If `start` = :cold, all links of the lattice are set to be the identity element.
If `start` = :hot, all links are randomized and normalized Sp2Element.
If `start` = :empty, the lattice contains initialized vectors of undef links.

`type` refers to what kind of representation to use for the Sp2 element in links:
- Sp2ElementA is the one used in Pepe's paper
- Sp2ElementB is the one used in Pepe's source code

"""
struct Lattice{D} <: AbstractArray{Vector{Link{D}}, D}
	lattice::Array{Vector{Link{D}}, D}
	
	function Lattice(dimensions::NTuple{D, Integer}; type::Type{<:Sp2Element} = Sp2ElementA, start::Symbol = :cold) where D
		if start ∉ [:cold, :hot, :empty]
			throw(ArgumentError("start keyword argument must be equal to :empty, :cold or :hot. Got start = :$start."))
		elseif D == 0
			throw(ArgumentError("At least one dimension must be provided to Lattice constructor."))
		elseif any(d -> d ≤ 0, dimensions)
			throw(ArgumentError("All dimensions provided to Lattice constructor must be strictly greater than zero."))
		end

		lattice = Array{Vector{Link}}(undef, dimensions)

		for index in CartesianIndices(lattice)
			lattice[index] = if start == :cold
				[Link(direction, index, type = type) for direction in 1:D]
			elseif start == :hot
				[Link(rand(type), direction, index) for direction in 1:D]
			else
				Vector{Link}(undef, D)
			end			
		end

		new{D}(lattice)
	end

	function Lattice(dimensions::Vararg{Integer, D}; type::Type{<:Sp2Element} = Sp2ElementA, start::Symbol = :cold) where D
		Lattice(Tuple(dimensions); type = type, start = start)
	end

	function Lattice(::Val{D}, length::Integer; type::Type{<:Sp2Element} = Sp2ElementA, start::Symbol = :cold) where D
		Lattice(ntuple(_ -> length, Val(D)); type = type, start = start)
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
Base.getindex(L::Lattice, i) = getindex(L.lattice, i)

# Multiple indices: every entry is `mod`ed by the dimensions of the lattice to wrap around
Base.getindex(L::Lattice, i...) = getindex(L.lattice, mod1.(i, size(L))...)

# Linear indexing, without `mod`
Base.setindex!(L::Lattice, v, i) = setindex!(L.lattice, v, i)

# Multiple indices: every entry is `mod`ed by the dimensions of the lattice to wrap around
Base.setindex!(L::Lattice, v, i...) = setindex!(L.lattice, v, mod1.(i, size(L))...)

Base.firstindex(L::Lattice) = firstindex(L.lattice)

Base.lastindex(L::Lattice) = lastindex(L.lattice)

Base.size(L::Lattice) = size(L.lattice)


#* ===== updatelattice! =====
"""
	updatelattice!(L::Lattice{D}, link::Link{D}) where D
Substitute the link already present in the lattice `L` with the provided `link`.
"""
updatelattice!(L::Lattice{D}, link::Link{D}) where D = L[link.position][link.direction] = link

"""
	updatelattice!(L::Lattice{D}, links::Vector{Link{D}}) where D
Substitute all links in `links` in the lattice `L`.
"""
function updatelattice!(L::Lattice{D}, links) where D
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
	nsites::Integer

	function EvenOddLattice(lattice::Lattice{D}, mod2_result::Int) where D
		if mod2_result ≠ 0 && mod2_result ≠ 1
			throw(ArgumentError("mod2_result in EvenOddLattice must be 0 or 1, got $mod2_result."))
		end
		nsites = count(p -> sum(p) % 2 == mod2_result, Tuple.(CartesianIndices(lattice)))
		new{D}(lattice, mod2_result, nsites)
	end
end

Base.iterate(L::EvenOddLattice) = iterate(L, Tuple(CartesianIndices(L.lattice)))

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

Base.length(L::EvenOddLattice) = L.nsites

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

"""
Iterator over a lattice links that points in the given `direction`.
"""
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

Base.iterate(L::EvenOddLinks) = iterate(L, Tuple(CartesianIndices(L.eolattice.lattice)))

function Base.iterate(L::EvenOddLinks, state)
	result = iterate(L.eolattice, state)
	if !isnothing(result)
		return (result[1][L.direction], result[2])
	end
	nothing
end

Base.length(L::EvenOddLinks) = L.eolattice.nsites

"""
Returns an iterator over all links pointing in the given `direction` that belong to even sites. See also `evensites`.
"""
evenlinks(L::Lattice, direction::Integer) = EvenOddLinks(evensites(L), direction)

"""
Returns an iterator over all links pointing in the given `direction` that belong to odd sites. See also `oddsites`.
"""
oddlinks(L::Lattice, direction::Integer) = EvenOddLinks(oddsites(L), direction)

function iterlinks(L::Lattice, direction::Integer, parity::Symbol) 
	if parity == :even
		evenlinks(L, direction)
	elseif parity == :odd
		oddlinks(L, direction)
	else
		throw(ArgumentError("iterlinks accepts either :even or :odd as third argument: got :$parity."))
	end
end

function everylink(L::Lattice{D}) where D
	Base.Iterators.flatten(
		[iterlinks(L, u, parity) for u in 1:D, parity in (:even, :odd)]
	)
end