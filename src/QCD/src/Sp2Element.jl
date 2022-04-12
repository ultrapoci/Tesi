export Sp2Element, asmatrix, normalizeSp2, randSp2

using LinearAlgebra

struct Sp2Element <: AbstractMatrix{ComplexF64}
	topleft::Matrix{ComplexF64}
	topright::Matrix{ComplexF64}

	function Sp2Element()
		new(Matrix{ComplexF64}(I, 2, 2), zeros(ComplexF64, 2, 2))
	end

	function Sp2Element(W::Matrix{ComplexF64}, X::Matrix{ComplexF64})
		if size(W) ≠ (2, 2) || size(X) ≠ (2, 2)
			throw(ArgumentError(
				"One or both of the given matrices don't respect the required dimension of (2,2):\n\
				\tsize(top left) = $(size(W)),\n\
				\tsize(top right) = $(size(X))."
			))
		end
		new(W, X)
	end

	function Sp2Element(W::Matrix{<:Number}, X::Matrix{<:Number})
		if size(W) ≠ (2, 2) || size(X) ≠ (2, 2)
			throw(ArgumentError(
				"One or both of the given matrices don't respect the required dimension of (2,2):\n\
				\tsize(top left) = $(size(W)),\n\
				\tsize(top right) = $(size(X))."
			))
		end
		new(ComplexF64.(W), ComplexF64.(X))
	end
end

"""
	normalizeSp2(S::Sp2Element)
Adjust the element of `S` to make sure that `S` belongs to SU(4), which means that `S` must be unitary and have determinant = 1.
"""
function normalizeSp2(S::Sp2Element)
	V₁ = [S.topleft[1, :]; S.topright[1, :]] # first row of S
	V₂ = [S.topleft[2, :]; S.topright[2, :]] # second row of S
	V₃ = [-conj(S.topright[1, :]); conj(S.topleft[1, :])] # third row of S

	# norm(V₁) = norm(V₃)
	N = norm(V₁)
	V₁ = V₁ / N
	V₃ = V₃ / N

	# NOTE: the scalar product ⋅ automatically takes the conjugate of the left term
	V₂ = V₂ - (V₁ ⋅ V₂)V₁ - (V₃ ⋅ V₂)V₃
	V₂ = V₂ / norm(V₂)

	W = [V₁[1] V₁[2]; V₂[1] V₂[2]]
	X = [V₁[3] V₁[4]; V₂[3] V₂[4]]

	Sp2Element(W, X)
end

function asmatrix(S::Sp2Element)
	[S.topleft S.topright; -conj(S.topright) conj(S.topleft)]
end

"""
Returns a random `Sp2Element`, already normalized such that it belongs to SU(4).
"""
function randSp2()
	normalizeSp2(Sp2Element(rand(ComplexF64, 2, 2), rand(ComplexF64, 2, 2)))
end

function Base.:*(S::Sp2Element, T::Sp2Element)
	R = asmatrix(S) * asmatrix(T)
	Sp2Element(R[1:2, 1:2], R[1:2, 3:4])
end

function Base.adjoint(S::Sp2Element)
	R = adjoint(asmatrix(S))
	Sp2Element(R[1:2, 1:2], R[1:2, 3:4])
end

function Base.getindex(S::Sp2Element, i...)
	getindex(asmatrix(S), i...)
end

function Base.setindex!(S::Sp2Element, v, i...)
	setindex!(asmatrix(S), v, i...)
end

function Base.firstindex(S::Sp2Element)
	firstindex(asmatrix(S))
end

function Base.lastindex(S::Sp2Element)
	lastindex(asmatrix(S))
end

function Base.size(S::Sp2Element)
	size(asmatrix(S))
end

function Base.zero(::Type{Sp2Element})
	z = zeros(ComplexF64, 2, 2)
	Sp2Element(z, z)
end

function LinearAlgebra.inv(S::Sp2Element)
	R = inv(asmatrix(S))
	Sp2Element(R[1:2, 1:2], R[1:2, 3:4])
end

function LinearAlgebra.tr(S::Sp2Element)
	tr(asmatrix(S))
end