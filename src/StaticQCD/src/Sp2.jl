export Sp2, asmatrix, normalizeSp2

struct Sp2 <: StaticMatrix{4, 4, ComplexF64}
	topleft::SMatrix{2, 2, ComplexF64, 4}
	topright::SMatrix{2, 2, ComplexF64, 4}

	Sp2(W::SMatrix{2, 2, ComplexF64, 4}, X::SMatrix{2, 2, ComplexF64, 4}) = new(W, X)

	Sp2() = Sp2(I[1:2, 1:2], zeros(2, 2))

	function Sp2(W::Matrix{<:Number}, X::Matrix{<:Number})
		Sp2(SMatrix{2, 2}(ComplexF64.(W)), SMatrix{2, 2}(ComplexF64.(X)))
	end
end

function asmatrix(S::Sp2)::SMatrix{4, 4, ComplexF64, 16}
	R = MMatrix{4, 4, ComplexF64}(undef)

	@inbounds R[1:2, 1:2] = S.topleft
	@inbounds R[1:2, 3:4] = S.topright

	@inbounds R[3, 1] =  conj(S.topright[2, 2])
	@inbounds R[3, 2] = -conj(S.topright[2, 1])
	@inbounds R[4, 1] = -conj(S.topright[1, 2])
	@inbounds R[4, 2] =  conj(S.topright[1, 1])

	@inbounds R[3, 3] =  conj(S.topleft[2, 2])
	@inbounds R[3, 4] = -conj(S.topleft[2, 1])
	@inbounds R[4, 3] = -conj(S.topleft[1, 2])
	@inbounds R[4, 4] =  conj(S.topleft[1, 1])

	SMatrix(R)
end

Base.getindex(S::Sp2, i::Int) = getindex(asmatrix(S), i)
Base.getindex(S::Sp2, i) = getindex(asmatrix(S), i)

function Base.:*(S::Sp2, T::Sp2)::Sp2
	R = asmatrix(S) * asmatrix(T)
	Sp2(R[1:2, 1:2], R[1:2, 3:4])
end

function Base.adjoint(S::Sp2)::Sp2
	R = adjoint(asmatrix(S))
	Sp2(R[1:2, 1:2], R[1:2, 3:4])
end

function LinearAlgebra.inv(S::Sp2)::Sp2
	R = inv(asmatrix(S))
	Sp2(R[1:2, 1:2], R[1:2, 3:4])
end

LinearAlgebra.tr(S::Sp2) = 2(real(S.topleft[1, 1]) + real(S.topleft[2, 2]))

Base.rand(::Type{Sp2}) = normalizeSp2(Sp2(rand(ComplexF64, 2, 2) .- complex(0.5, 0.5), rand(ComplexF64, 2, 2) .- complex(0.5, 0.5)))

function normalizeSp2(S::Sp2)::Sp2
	s = asmatrix(S)
	@inbounds V₁ = s[1, :] # first row of S
	@inbounds V₂ = s[2, :] # second row of S
	@inbounds V₃ = s[4, :] # fourth row of S

	# norm(V₁) = norm(V₃)
	N = norm(V₁)
	V₁ = V₁ / N
	V₃ = V₃ / N

	# NOTE: the scalar product ⋅ automatically takes the conjugate of the left term
	V₂ = V₂ - (V₁ ⋅ V₂)V₁ - (V₃ ⋅ V₂)V₃
	V₂ = V₂ / norm(V₂)

	@inbounds W = [V₁[1] V₁[2]; V₂[1] V₂[2]]
	@inbounds X = [V₁[3] V₁[4]; V₂[3] V₂[4]]

	Sp2(W, X)
end