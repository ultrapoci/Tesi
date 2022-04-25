export SU2Element, asmatrix, normalizeSU2, normalizeSU2det

struct SU2Element <: AbstractMatrix{ComplexF64}
	t₁::ComplexF64
	t₂::ComplexF64

	function SU2Element(t₁::Number, t₂::Number)
		SU2Element(ComplexF64(t₁), ComplexF64(t₂))
	end

	function SU2Element(t₁::ComplexF64, t₂::ComplexF64)
		new(t₁, t₂)
	end
end

asmatrix(S::SU2Element) = [S.t₁ S.t₂; -conj(S.t₂) conj(S.t₁)]

function normalizeSU2(S::SU2Element)
	sqrtΔ = √det(S)
	t₁ = S.t₁ / sqrtΔ
	t₂ = S.t₂ / sqrtΔ
	SU2Element(t₁, t₂)
end

function normalizeSU2det(S::SU2Element)
	sqrtΔ = √det(S)
	t₁ = S.t₁ / sqrtΔ
	t₂ = S.t₂ / sqrtΔ
	SU2Element(t₁, t₂), sqrtΔ
end

function Base.:*(S::SU2Element, T::SU2Element)
	R = asmatrix(S) * asmatrix(T)
	SU2Element(R[1, 1], R[1, 2])
end

function Base.:+(S::SU2Element, T::SU2Element)
	R = asmatrix(S) + asmatrix(T)
	SU2Element(R[1, 1], R[1, 2])
end

function Base.adjoint(S::SU2Element)
	R = adjoint(asmatrix(S))
	SU2Element(R[1, 1], R[1, 2])
end

Base.getindex(S::SU2Element, i...) = getindex(asmatrix(S), i...)

Base.setindex!(S::SU2Element, v, i...) = setindex!(asmatrix(S), v, i...)

Base.firstindex(S::SU2Element) = firstindex(asmatrix(S))

Base.lastindex(S::SU2Element) = lastindex(asmatrix(S))

Base.size(S::SU2Element) = size(asmatrix(S))

Base.zero(::Type{SU2Element}) = SU2Element(zero(ComplexF64), zero(ComplexF64))

LinearAlgebra.det(S::SU2Element) = abs2(S.t₁) + abs2(S.t₂)

function LinearAlgebra.inv(S::SU2Element)
	R = inv(asmatrix(S))
	SU2Element(R[1, 1], R[1, 2])
end

LinearAlgebra.tr(S::SU2Element) = 2 * real(S.t₁)
