export Sp2Element, asmatrix

# TODO check if it is also an element of SU(4)
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

end

function asmatrix(S::Sp2Element)
	[S.topleft S.topright; -conj(S.topright) conj(S.topleft)]
end

function Base.:*(S::Sp2Element, T::Sp2Element)
	asmatrix(S) * asmatrix(T)
end

function Base.adjoint(S::Sp2Element)
	adjoint(asmatrix(S))
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
end