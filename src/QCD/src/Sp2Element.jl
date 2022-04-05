export Sp2Element, asmatrix

# TODO check if it is also an element of SU(4)
struct Sp2Element <: AbstractArray{Number, 2}
	topleft::Matrix{<:Number}
	topright::Matrix{<:Number}

	function Sp2Element()
		new(Matrix{Complex{Int}}(I, 2, 2), zeros(Complex{Int}, 2, 2))
	end

	function Sp2Element(W::Matrix{<:Number}, X::Matrix{<:Number})
		if size(W) ≠ (2, 2) || size(X) ≠ (2, 2)
			throw(
				ArgumentError(
					"One or both of the given matrices don't respect the required dimension of (2,2):\n\tsize(W) = $(size(W)),\n\tsize(X) = $(size(X))."
				)
			)	
		end
		new(W, X)
	end
end

function asmatrix(S::Sp2Element)
	[S.topleft S.topright; -conj(S.topright) conj(S.topleft)]
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
