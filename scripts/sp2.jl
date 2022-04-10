using LinearAlgebra, QCD

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