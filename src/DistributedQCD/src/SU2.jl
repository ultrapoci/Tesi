struct SU2 <: StaticMatrix{2, 2, ComplexF64}
	t₁::ComplexF64
	t₂::ComplexF64

	SU2(t₁::Number, t₂::Number) = SU(ComplexF64(t₁), ComplexF64(t₂))
	SU2(t₁::ComplexF64, t₂::ComplexF64) = new(t₁, t₂)
end

asmatrix(S::SU2) = SMatrix{2, 2, ComplexF64}([S.t₁ S.t₂; -conj(S.t₂) conj(S.t₁)])

Base.getindex(S::SU2, i::Int) = getindex(asmatrix(S), i)
Base.getindex(S::SU2, i) = getindex(asmatrix(S), i)

Base.:*(S::SU2, T::SU2) = SU2(S.t₁*T.t₁ - S.t₂*conj(T.t₂), S.t₁*T.t₂ + S.t₂*conj(T.t₁))
Base.adjoint(S::SU2) = SU2(conj(S.t₁), -S.t₂)

LinearAlgebra.det(S::SU2) = abs2(S.t₁) + abs2(S.t₂)
LinearAlgebra.tr(S::SU2) = 2real(S.t₁)
function LinearAlgebra.inv(S::SU2)
	Δ = det(S)
	SU2(conj(S.t₁) / Δ, -S.t₂ / Δ)
end

function normalizeSU2(S::SU2)
	sqrtΔ = √det(S)
	t₁ = S.t₁ / sqrtΔ
	t₂ = S.t₂ / sqrtΔ
	SU2(t₁, t₂)
end

function normalizeSU2det(S::SU2)
	sqrtΔ = √det(S)
	t₁ = S.t₁ / sqrtΔ
	t₂ = S.t₂ / sqrtΔ
	SU2(t₁, t₂), sqrtΔ
end
