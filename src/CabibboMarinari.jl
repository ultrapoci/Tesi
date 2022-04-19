using LinearAlgebra, QCD
using Distributions: Uniform

function iterlinks(L::Lattice, direction::Integer, parity::Symbol) 
	if parity == :even
		evenlinks(L, direction)
	elseif parity == :odd
		oddlinks(L, direction)
	else
		throw(ArgumentError("iterlinks accepts either :even or :odd as third argument: got :$parity."))
	end
end

function sumstaples(lattice::Lattice{D}, link::Link{D}) where D
	total = zeros(ComplexF64, 4, 4)
	u = link.direction
	for v in mod1.(u+1:u+D-1, D) # generate all D dimensions except u
		s₊ = staple(lattice, link, v) |> asmatrix
		s₋ = staple(lattice, link, -v)  |> asmatrix
		total += s₊ + s₋
	end
	total
end

"""
	subrepresentations(::Type{Sp2ElementA})
	subrepresentations(::Type{Sp2ElementB})
Returns a vector of tuples, each containing two functions: 
- the first function takes an Sp2Element and returns a tuple, where the first element is the extracted SU2Element (already normalized) according to \
one of the four decompositions of the ``Sp2`` algebra, while the second element is the square root of the determinant calculated before the normalization
- the second function does the inverse: takes an SU2Element and builds an Sp2Element according the decomposition considered. 

Each pair represent one algebra decomposition. 

The input type determines which form of the Sp2 matrix to use when building decompositions.
"""
function subrepresentations(::Type{Sp2ElementA})
	repr = Tuple{Function, Function}[]

	# SU(2) + U(1) + U(1)
	push!(
		repr, 
		(
			x -> SU2Element(x[1, 1], x[1, 3]) |> normalizeSU2det,
			(x::SU2Element) -> Sp2Element(Sp2ElementA, [x.t₁ 0; 0 1], [x.t₂ 0; 0 0])
		)
	)

	# SU(2) + U(1) + U(1)
	push!(
		repr, 
		(
			x -> SU2Element(x[2, 2], x[2, 4]) |> normalizeSU2det,
			(x::SU2Element) -> Sp2Element(Sp2ElementA, [1 0; 0 x.t₁], [0 0; 0 x.t₂])
		)
	)

	# SU(2) + SU(2)
	push!(
		repr, 
		(
			x -> SU2Element(x[1, 1] + x[2, 2], x[1, 4] + x[2, 3]) |> normalizeSU2det,
			(x::SU2Element) -> Sp2Element(Sp2ElementA, [x.t₁ 0; 0 x.t₁], [0 x.t₂; x.t₂ 0])
		)
	)

	# SU(2) + SU(2)
	push!(
		repr, 
		(
			x -> SU2Element(x[1, 1] + x[4, 4], x[1, 2] - x[4, 3]) |> normalizeSU2det,
			(x::SU2Element) -> Sp2Element(Sp2ElementA, [x.t₁ x.t₂; -conj(x.t₂) conj(x.t₁)], [0 0; 0 0])
		)
	)

	repr
end

function subrepresentations(::Type{Sp2ElementB})
	repr = Tuple{Function, Function}[]

	# SU(2) + U(1) + U(1)
	push!(
		repr, 
		(
			x -> SU2Element(x[1, 1], x[1, 4]) |> normalizeSU2det,
			(x::SU2Element) -> Sp2Element(Sp2ElementB, [x.t₁ 0; 0 1], [0 x.t₂; 0 0])
		)
	)

	# SU(2) + U(1) + U(1)
	push!(
		repr, 
		(
			x -> SU2Element(x[2, 2], x[2, 3]) |> normalizeSU2det,
			(x::SU2Element) -> Sp2Element(Sp2ElementB, [1 0; 0 x.t₁], [0 0; x.t₂ 0])
		)
	)

	# SU(2) + SU(2)
	push!(
		repr, 
		(
			x -> SU2Element(x[1, 1] + x[2, 2], x[1, 3] - x[2, 4]) |> normalizeSU2det,
			(x::SU2Element) -> Sp2Element(Sp2ElementB, [x.t₁ 0; 0 x.t₁], [x.t₂ 0; 0 -x.t₂])
		)
	)

	# SU(2) + SU(2)
	push!(
		repr, 
		(
			x -> SU2Element(x[1, 1] + x[3, 3], x[1, 2] + x[3, 4]) |> normalizeSU2det,
			(x::SU2Element) -> Sp2Element(Sp2ElementB, [x.t₁ x.t₂; -conj(x.t₂) conj(x.t₁)], [0 0; 0 0])
		)
	)

	repr
end

"""
	generate_a0(k::Real, β::Real)
Generates a real number ``a₀`` according to the distribution P(a₀) = √(1 - a₀^2) exp(a₀ β k).		
"""
function generate_a0(k::Real, β::Real)
	reject = true
	a₀ = 0.0
	while reject
		x = rand(Uniform(exp(-2*β*k), 1.0)) #? are range extremes a problem?
		a₀ = 1 + log(x) / (β*k)
		reject = 1 - √(1 - a₀^2) > rand(Uniform(0.0, 1.0))
	end	
	a₀
end

function randomSU2(k::Real, β::Real)
	a₀ = generate_a0(k, β)
	ϕ = rand(Uniform(0.0, 2π))
	θ = acos(rand(Uniform(-1.0, 1.0)))

	r = √(1 - a₀^2)

	a₁ = r * sin(θ) * cos(ϕ)
	a₂ = r * sin(θ) * sin(ϕ)
	a₃ = r * cos(θ)

	SU2Element(complex(a₀, a₃), complex(a₂, a₁))
end

function overrelaxation(lattice::Lattice, link::Link)
	U = link.s
	T = typeof(U)
	R = sumstaples(lattice, link)
	for (to_su2, to_sp2) in subrepresentations(T)
		u, = to_su2(U * R) # don't need the square root of the determinant
		a = u^-2
		U = to_sp2(a) * U
	end
	Link(U, link.direction, link.position)
end

function heatbath(lattice::Lattice, link::Link, β::Real)
	U = link.s
	T = typeof(U)
	R = sumstaples(lattice, link)
	for (to_su2, to_sp2) in subrepresentations(T)
		u, k = to_su2(U * R)
		a = randomSU2(k, β) * u^-1
		U = to_sp2(a) * U
	end
	Link(U, link.direction, link.position)
end

function lattice_overrelaxation!(lattice::Lattice{D}, n::Integer) where D
	for _ in 1:n, parity in [:even, :odd], u in 1:D
		newlinks = Link{D}[]
		for link in iterlinks(lattice, u, parity)
			newlink = overrelaxation(lattice, link)
			push!(newlinks, newlink)
		end
		updatelattice!(lattice, newlinks)
	end
end

function lattice_heatbath!(lattice::Lattice{D}, β::Real) where D
	for parity in [:even, :odd], u in 1:D
		newlinks = Link{D}[]
		for link in iterlinks(lattice, u, parity)
			newlink = heatbath(lattice, link, β)
			push!(newlinks, newlink)
		end
		updatelattice!(lattice, newlinks)
	end
end

function lattice_normalization!(lattice::Lattice{D}) where D
	for x in CartesianIndices(lattice), u in 1:D
		s = normalizeSp2(linkelement(lattice, u, x))
		updatelattice!(lattice, Link(s, u, x))
	end
end

function averageplaquette(lattice::Lattice{D}) where D
	s = 0.0
	for site in lattice, link in site[begin:end-1], v in link.direction+1:D
		s += tr(plaquette(lattice, link, v))
	end
	
	V = length(lattice) # lattice's volume
	np = D * (D - 1) ÷ 2 # number of plaquettes per site
	
	# the 4 term is due to 2N for N=2
	s / (4 * np * V) 
end

function action(lattice::Lattice{D}) where D
	s = 0.0
	for site in lattice, link in site[begin:end-1], v in link.direction+1:D
		s += 1 - 0.5 * tr(plaquette(lattice, link, v))
	end
	s
end

function polyakovloop(lattice::Lattice{D}, x::CartesianIndex{Dm1}) where {D, Dm1}
	Dm1 == D-1 || throw(TypeError(:polyakovloop, CartesianIndex{D-1}, CartesianIndex{Dm1}))
	U = linkelement(lattice, 1, CartesianIndex(1, x))
	Nₜ = first(size(lattice))
	for i in 2:Nₜ
		U = U * linkelement(lattice, 1, CartesianIndex(i, x))
	end
	tr(U)
end

polyakovloop(lattice::Lattice, t::Tuple{Vararg{Integer}}) = polyakovloop(lattice, CartesianIndex(t))
polyakovloop(lattice::Lattice, t::Vararg{Integer}) = polyakovloop(lattice, CartesianIndex(t))