export Sp2Element, Sp2ElementA, Sp2ElementB, asmatrix, normalizeSp2, J

abstract type Sp2Element <: AbstractMatrix{ComplexF64} end

struct Sp2ElementA <: Sp2Element
	topleft::Matrix{ComplexF64}
	topright::Matrix{ComplexF64}

	function Sp2ElementA()
		new(Matrix{ComplexF64}(I, 2, 2), zeros(ComplexF64, 2, 2))
	end

	function Sp2ElementA(W::Matrix{ComplexF64}, X::Matrix{ComplexF64})
		if size(W) ≠ (2, 2) || size(X) ≠ (2, 2)
			throw(ArgumentError(
				"One or both of the given matrices don't respect the required dimension of (2,2):\n\
				\tsize(top left) = $(size(W)),\n\
				\tsize(top right) = $(size(X))."
			))
		end
		new(W, X)
	end

	function Sp2ElementA(W::Matrix{<:Number}, X::Matrix{<:Number})
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

# from Pepe's source code
struct Sp2ElementB <: Sp2Element
	topleft::Matrix{ComplexF64}
	topright::Matrix{ComplexF64}

	function Sp2ElementB()
		new(Matrix{ComplexF64}(I, 2, 2), zeros(ComplexF64, 2, 2))
	end

	function Sp2ElementB(W::Matrix{ComplexF64}, X::Matrix{ComplexF64})
		if size(W) ≠ (2, 2) || size(X) ≠ (2, 2)
			throw(ArgumentError(
				"One or both of the given matrices don't respect the required dimension of (2,2):\n\
				\tsize(top left) = $(size(W)),\n\
				\tsize(top right) = $(size(X))."
			))
		end
		new(W, X)
	end

	function Sp2ElementB(W::Matrix{<:Number}, X::Matrix{<:Number})
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

Sp2Element(::Type{Sp2ElementA}) = Sp2ElementA()
Sp2Element(::Type{Sp2ElementB}) = Sp2ElementB()
Sp2Element(::Type{Sp2ElementA}, W, X) = Sp2ElementA(W, X)
Sp2Element(::Type{Sp2ElementB}, W, X) = Sp2ElementB(W, X)

J(::Type{Sp2ElementA}) = [0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0]
J(::Type{Sp2ElementB}) = [0 0 0 1; 0 0 -1 0; 0 1 0 0; -1 0 0 0]

asmatrix(S::Sp2ElementA) = [S.topleft S.topright; -conj(S.topright) conj(S.topleft)]

function asmatrix(S::Sp2ElementB)
	X = conj.([S.topright[2, 2] -S.topright[2, 1]; -S.topright[1, 2] S.topright[1, 1]])
	W = conj.([S.topleft[2, 2] -S.topleft[2, 1]; -S.topleft[1, 2] S.topleft[1, 1]])
	[S.topleft S.topright; X W]
end

function Base.rand(T::Type{<:Sp2Element})
	normalizeSp2(
		Sp2Element(
			T,
			rand(ComplexF64, 2, 2) .- complex(0.5, 0.5), 
			rand(ComplexF64, 2, 2) .- complex(0.5, 0.5)
		)
	)
end

"""
	normalizeSp2(S::Sp2Element)
Adjust the element of `S` to make sure that `S` belongs to SU(4), which means that `S` must be unitary and have determinant = 1.
"""
function normalizeSp2(S::Sp2ElementA)
	s = asmatrix(S)
	V₁ = s[1, :] # first row of S
	V₂ = s[2, :] # second row of S
	V₃ = s[3, :] # third row of S

	# norm(V₁) = norm(V₃)
	N = norm(V₁)
	V₁ = V₁ / N
	V₃ = V₃ / N

	# NOTE: the scalar product ⋅ automatically takes the conjugate of the left term
	V₂ = V₂ - (V₁ ⋅ V₂)V₁ - (V₃ ⋅ V₂)V₃
	V₂ = V₂ / norm(V₂)

	W = [V₁[1] V₁[2]; V₂[1] V₂[2]]
	X = [V₁[3] V₁[4]; V₂[3] V₂[4]]

	Sp2ElementA(W, X)
end

function normalizeSp2(S::Sp2ElementB)
	s = asmatrix(S)
	V₁ = s[1, :] # first row of S
	V₂ = s[2, :] # second row of S
	V₃ = s[4, :] # fourth row of S

	# norm(V₁) = norm(V₃)
	N = norm(V₁)
	V₁ = V₁ / N
	V₃ = V₃ / N

	# NOTE: the scalar product ⋅ automatically takes the conjugate of the left term
	V₂ = V₂ - (V₁ ⋅ V₂)V₁ - (V₃ ⋅ V₂)V₃
	V₂ = V₂ / norm(V₂)

	W = [V₁[1] V₁[2]; V₂[1] V₂[2]]
	X = [V₁[3] V₁[4]; V₂[3] V₂[4]]

	Sp2ElementB(W, X)
end

function Base.:*(S::Sp2Type, T::Sp2Type) where Sp2Type <: Sp2Element
	R = asmatrix(S) * asmatrix(T)
	Sp2Element(Sp2Type, R[1:2, 1:2], R[1:2, 3:4])
end

Base.:*(S::Sp2Element, T::Matrix{<:Number}) = asmatrix(S) * T

Base.:*(S::Matrix{<:Number}, T::Sp2Element) = S * asmatrix(T)

function Base.adjoint(S::Sp2Element)
	R = adjoint(asmatrix(S))
	Sp2Element(typeof(S), R[1:2, 1:2], R[1:2, 3:4])
end

Base.getindex(S::Sp2Element, i...) = getindex(asmatrix(S), i...)

Base.setindex!(S::Sp2Element, v, i...) = setindex!(asmatrix(S), v, i...)

Base.firstindex(S::Sp2Element) = firstindex(asmatrix(S))

Base.lastindex(S::Sp2Element) = lastindex(asmatrix(S))

Base.size(::Sp2Element) = (4, 4)

function Base.zero(T::Type{<:Sp2Element})
	z = zeros(ComplexF64, 2, 2)
	Sp2Element(T, z, z)
end

function LinearAlgebra.inv(S::Sp2Element)
	R = inv(asmatrix(S))
	Sp2Element(typeof(S), R[1:2, 1:2], R[1:2, 3:4])
end

LinearAlgebra.tr(S::Sp2Element) = 2(real(S.topleft[1, 1]) + real(S.topleft[2, 2]))