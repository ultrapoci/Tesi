using LinearAlgebra, QCD, Distributions

function staplesum(lattice::Lattice{D}, link::Link{D}) where D
	total = zero(Sp2Element)
	u = link.direction
	for v in mod1.(u+1:u+D-1, D) # generate all D dimensions except u
		s₊ = getstaple(lattice, link, v)
		s₋ = getstaple(lattice, link, -v)
		total += s₊ + s₋
	end
	total
end

function subrepresentations(s::Sp2Element)
	S = asmatrix(s)
	repr = Tuple{Matrix{Complex}, Real}[]

	# SU(2) + U(1) + U(1)
	t₁ = S[1, 1]
	t₂ = S[1, 3]
	push!(repr, SU2Element(t₁, t₂) |> normalizeSU2det)

	# SU(2) + U(1) + U(1)
	t₁ = S[2, 2]
	t₂ = S[2, 4]
	push!(repr, SU2Element(t₁, t₂) |> normalizeSU2det)

	# SU(2) + SU(2)
	t₁ = S[1, 1] + S[2, 2]
	t₂ = S[1, 4] + S[2, 3]
	push!(repr, SU2Element(t₁, t₂) |> normalizeSU2det)

	# SU(2) + SU(2)
	t₁ = S[1, 1] + S[4, 4]
	t₂ = S[1, 2] - S[4, 3]
	push!(repr, SU2Element(t₁, t₂) |> normalizeSU2det)

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