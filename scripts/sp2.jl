using LinearAlgebra, QCD, Distributions

function sumstaples(lattice::Lattice{D}, link::Link{D}) where D
	total = zero(Sp2Element)
	u = link.direction
	for v in mod1.(u+1:u+D-1, D) # generate all D dimensions except u
		s₊ = staple(lattice, link, v)
		s₋ = staple(lattice, link, -v)
		total += s₊ + s₋
	end
	total
end

"""
	subrepresentations(s::Sp2Element)
Extracts the four ``SU(2)`` matrices embedded in an ``Sp(2)`` element and returns a tuple containing
the SU2Element already normalized, the square root of the determinant calculated before the normalization
and an anonymous function that takes the two complex numbers in SU2Element and builds an Sp2Element.
"""
function subrepresentations(s::Sp2Element)
	S = asmatrix(s)
	repr = Tuple{SU2Element, Real, Function}[]

	# SU(2) + U(1) + U(1)
	t₁ = S[1, 1]
	t₂ = S[1, 3]
	(M, sqrtΔ) = normalizeSU2det(SU2Element(t₁, t₂))
	push!(
		repr, 
		(M, sqrtΔ, (x, y) -> Sp2Element([x 0; 0 1], [y 0; 0 0]))
	)

	# SU(2) + U(1) + U(1)
	t₁ = S[2, 2]
	t₂ = S[2, 4]
	(M, sqrtΔ) = normalizeSU2det(SU2Element(t₁, t₂))
	push!(
		repr, 
		(M, sqrtΔ, (x, y) -> Sp2Element([1 0; 0 x], [0 0; 0 y]))
	)

	# SU(2) + SU(2)
	t₁ = S[1, 1] + S[2, 2]
	t₂ = S[1, 4] + S[2, 3]
	(M, sqrtΔ) = normalizeSU2det(SU2Element(t₁, t₂))
	push!(
		repr, 
		(M, sqrtΔ, (x, y) -> Sp2Element([x 0; 0 x], [0 y; y 0]))
	)

	# SU(2) + SU(2)
	t₁ = S[1, 1] + S[4, 4]
	t₂ = S[1, 2] - S[4, 3]
	(M, sqrtΔ) = normalizeSU2det(SU2Element(t₁, t₂))
	push!(
		repr, 
		(M, sqrtΔ, (x, y) -> Sp2Element([x y; -conj(y) conj(x)], [0 0; 0 0]))
	)

	repr
end

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
		U = f(a.t₁, a.t₂) * U
	end
	U
end

function heatbath(lattice::Lattice, link::Link, β::Real)
	U = link.s
	R = sumstaples(lattice, link)
	for (u, k, f) in subrepresentations(U * R)
		a = randomSU2(k, β) * u^-1
		U = f(a.t₁, a.t₂) * U
	end
	U
end

function averageplaquette(lattice::Lattice{D}) where D
	s = 0.0
	for site in lattice, link in site[begin:end-1], v in link.direction+1:D
		s += tr(plaquette(lattice, link, v))
	end
	
	V = length(lattice)

	# the 4 term is due to 2N for N=2
	# sum(1:D-1) is the number of plaquettes per site
	s / (4 * sum(1:D-1) * V) 
end

function lattice_overrelaxation(lattice::Lattice{D}, n::Integer) where D
	dims = size(lattice)
	prevlattice = lattice
	for _ in 1:n
		newlattice = Lattice(dims, start = LatticeStart.Empty)
		#for itersites in [evensites(prevlattice), oddsites(prevlattice)], site in itersites, u in 1:D
		for iterlinks in [evenlink, oddlinks], u in 1:D, link in iterlinks(prevlattice, u)
			direction = link.direction
			position = link.position
			U = overrelaxation(lattice, link)
			newlattice[position][direction] = Link(U, direction, position)
		end
		prevlattice = newlattice
	end
	prevlattice
end