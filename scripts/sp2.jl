using LinearAlgebra, QCD

function SU2matrix(t₁::Complex, t₂::Complex)
	Δ = √(abs2(t₁) + abs2(t₂))
	m = [t₁ t₂; -conj(t₂) conj(t₁)]
	(m ./ Δ, Δ)
end

function subrepresentations(s::Sp2Element)
	S = asmatrix(s)
	repr = Tuple{Matrix{Complex}, Real}[]

	# SU(2) + U(1) + U(1)
	t₁ = S[1, 1]
	t₂ = S[1, 3]
	push!(repr, SU2matrix(t₁, t₂))

	# SU(2) + U(1) + U(1)
	t₁ = S[2, 2]
	t₂ = S[2, 4]
	push!(repr, SU2matrix(t₁, t₂))

	# SU(2) + SU(2)
	t₁ = S[1, 1] + S[2, 2]
	t₂ = S[1, 4] + S[2, 3]
	push!(repr, SU2matrix(t₁, t₂))

	# SU(2) + SU(2)
	t₁ = S[1, 1] + S[4, 4]
	t₂ = S[1, 2] - S[4, 3]
	push!(repr, SU2matrix(t₁, t₂))

	repr
end