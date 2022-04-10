export SU2Element, asmatrix, normalizeSU2, normalizeSU2det

using LinearAlgebra

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

function asmatrix(S::SU2Element)
	[S.t₁ S.t₂; -conj(S.t₂) conj(S.t₁)]
end

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

function Base.adjoint(S::SU2Element)
	R = adjoint(asmatrix(S))
	SU2Element(R[1, 1], R[1, 2])
end

function Base.getindex(S::SU2Element, i...)
	getindex(asmatrix(S), i...)
end

function Base.setindex!(S::SU2Element, v, i...)
	setindex!(asmatrix(S), v, i...)
end

function Base.firstindex(S::SU2Element)
	firstindex(asmatrix(S))
end

function Base.lastindex(S::SU2Element)
	lastindex(asmatrix(S))
end

function Base.size(S::SU2Element)
	size(asmatrix(S))
end

function LinearAlgebra.det(S::SU2Element)
	abs2(S.t₁) + abs2(S.t₂)
end