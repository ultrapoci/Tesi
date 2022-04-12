using LinearAlgebra, QCD, Distributions

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
	subrepresentations(S::Matrix{<:Number})
	subrepresentations(S::Sp2Element)
Extracts the four ``SU(2)`` matrices embedded in a 4x4 matrix and returns a tuple containing
the SU2Element already normalized, the square root of the determinant calculated before the normalization
and an anonymous function that takes the two complex numbers in SU2Element and builds an Sp2Element.
"""
function subrepresentations(S::Matrix{<:Number})
	repr = Tuple{SU2Element, Real, Function}[]

	# SU(2) + U(1) + U(1)
	t₁ = S[1, 1]
	t₂ = S[1, 3]
	(M, sqrtΔ) = normalizeSU2det(SU2Element(t₁, t₂))
	push!(
		repr, 
		(M, sqrtΔ, (x::SU2Element) -> Sp2Element([x.t₁ 0; 0 1], [x.t₂ 0; 0 0]))
	)

	# SU(2) + U(1) + U(1)
	t₁ = S[2, 2]
	t₂ = S[2, 4]
	(M, sqrtΔ) = normalizeSU2det(SU2Element(t₁, t₂))
	push!(
		repr, 
		(M, sqrtΔ, (x::SU2Element) -> Sp2Element([1 0; 0 x.t₁], [0 0; 0 x.t₂]))
	)

	# SU(2) + SU(2)
	t₁ = S[1, 1] + S[2, 2]
	t₂ = S[1, 4] + S[2, 3]
	(M, sqrtΔ) = normalizeSU2det(SU2Element(t₁, t₂))
	push!(
		repr, 
		(M, sqrtΔ, (x::SU2Element) -> Sp2Element([x.t₁ 0; 0 x.t₁], [0 x.t₂; x.t₂ 0]))
	)

	# SU(2) + SU(2)
	t₁ = S[1, 1] + S[4, 4]
	t₂ = S[1, 2] - S[4, 3]
	(M, sqrtΔ) = normalizeSU2det(SU2Element(t₁, t₂))
	push!(
		repr, 
		(M, sqrtΔ, (x::SU2Element) -> Sp2Element([x.t₁ x.t₂; -conj(x.t₂) conj(x.t₁)], [0 0; 0 0]))
	)

	repr
end

function subrepresentations(S::Sp2Element)
	subrepresentations(asmatrix(S))
end

"""
	generate_a0(k::Real, β::Real)
Generates a real number ``a₀`` according to the distribution √(1 - a₀^2) exp(a₀ β k).		
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
	R = sumstaples(lattice, link)
	for (u, _, f) in subrepresentations(U * R)
		a = u^-2
		U = f(a) * U
	end
	U
end

function heatbath(lattice::Lattice, link::Link, β::Real)
	U = link.s
	R = sumstaples(lattice, link)
	for (u, k, f) in subrepresentations(U * R)
		a = randomSU2(k, β) * u^-1
		U = f(a) * U
	end
	U
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

function lattice_overrelaxation!(lattice::Lattice{D}, n::Integer) where D
	for _ in 1:n, parity in [:even, :odd], u in 1:D
		newlinks = Link{D}[]
		for link in iterlinks(lattice, u, parity)
			U = overrelaxation(lattice, link)
			push!(newlinks, Link(U, link.direction, link.position))
		end
		updatelattice!(lattice, newlinks)
	end
end

function lattice_heatbath!(lattice::Lattice{D}, β::Real) where D
	for parity in [:even, :odd], u in 1:D
		newlinks = Link{D}[]
		for link in iterlinks(lattice, u, parity)
			U = heatbath(lattice, link, β)
			push!(newlinks, Link(U, link.direction, link.position))
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