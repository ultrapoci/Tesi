using LinearAlgebra, QCD, Distributions

const λ = [
	[1  0  0  0
	 0  0  0  0
	 0  0 -1  0
	 0  0  0  0], 

	[0  0  0  0
	 0  1  0  0
	 0  0  0  0
	 0  0  0 -1], 

	[0  1  0  0
	 1  0  0  0
	 0  0  0 -1
	 0  0 -1  0], 

	[0 im  0  0
	-im 0  0  0
	 0  0  0 im
	 0  0 -im 0],

	[0  0  1  0
	 0  0  0  0
	 1  0  0  0
	 0  0  0  0], 

	[0  0 im  0
	 0  0  0  0
	-im 0  0  0
	 0  0  0  0], 

	[0  0  0  0
	 0  0  0  1
	 0  0  0  0
	 0  1  0  0], 

	[0  0  0  0
	 0  0  0 im
	 0  0  0  0
	 0 -im 0  0], 

	[0  0  0  1
	 0  0  1  0
	 0  1  0  0
	 1  0  0  0], 

	[0  0  0 im
	 0  0 im  0
	 0 -im 0  0
	-im 0  0  0], 
]

const id = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
const J = [0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0]

const V = [
	[1, 0, 0, 0],
	[0, 1, 0, 0],
	[0, 0, 1, 0],
	[0, 0, 0, 1]
]

##

function randomSp2()
	v = im .* λ
	b = rand(Uniform(-1.0, 1.0), length(v))
	a = rand(Uniform(-1.0, 1.0), 2)
	diagm([a[1], a[2], a[1], a[2]]) + sum(b .* v)
end

function normalizematrixSp2(S::Matrix{<:Number})
	V₁ = S[1, :] # first row of S
	V₂ = S[2, :] # second row of S
	V₃ = S[3, :] # third row of S

	# norm(V₁) = norm(V₃)
	N = norm(V₁)
	V₁ = V₁ / N
	V₃ = V₃ / N

	# NOTE: the scalar product ⋅ automatically takes the conjugate of the left term
	V₂ = V₂ - (V₁ ⋅ V₂)V₁ - (V₃ ⋅ V₂)V₃
	V₂ = V₂ / norm(V₂)

	W = [V₁[1] V₁[2]; V₂[1] V₂[2]]
	X = [V₁[3] V₁[4]; V₂[3] V₂[4]]

	[W X; -conj(X) conj(W)]
end