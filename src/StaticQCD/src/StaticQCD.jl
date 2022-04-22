module StaticQCD

export newlattice, Site, LocalLattice, Lattice, initprocs, getlink, addtuple, evenmask

using Distributed, DistributedArrays, StaticArrays
include("StaticSp2.jl")

const Site{D} = MVector{D, Sp2}
const LocalLattice{D} = Array{Site{D}, D}
const Lattice{D} = DArray{Site{D}, D, LocalLattice{D}}

initprocs(n; kwargs...) = addprocs(n, exeflags="--project", kwargs...)

function newlattice(dims::NTuple{D, Int}; start = :cold) where D
	if start == :cold
		distribute([MVector{D}([Sp2() for _ in 1:D]) for _ in CartesianIndices(dims)])
	elseif start == :hot
		distribute([MVector{D}([rand(Sp2) for _ in 1:D]) for _ in CartesianIndices(dims)])
	else
		throw(ArgumentError("start must be :hot or :cold, got :$start."))
	end
end
newlattice(dims::Vararg{Int, D}; kwargs...) where D = newlattice(Tuple(dims); kwargs...)

addtuple(n::Int, d::Int, p::Tuple{Vararg{Int}}) = tuple(p[begin:d-1]..., p[d]+n, p[d+1:end]...)
addtuple(n::Int, d::Int, p::Vararg{Int}) = addtuple(n, d, Tuple(p))

function getlink(L::Lattice{D}, d::Int, p::NTuple{D, Int}) where D
	if d ∈ 1:D
		x = mod1.(p, size(L))
		L[x...][d]
	elseif d ∈ -D:-1
		x = mod1.(addtuple(-1, -d, p), size(L))
		L[x...][-d]'
	else
		throw(ArgumentError("Link direction must be in range [1, $D] or [-$D, -1]. Got d = $d"))
	end
end
getlink(L::Lattice{D}, d::Int, p::Vararg{Int, D}) where D = getlink(L, d, Tuple(p))
getlink(L::Lattice{D}, d::Int, p::CartesianIndex{D}) where D = getlink(L, d, Tuple(p))

function getlink(L::LocalLattice{D}, d::Int, p::NTuple{D, Int}) where D
	if d ∈ 1:D
		x = mod1.(p, size(L))
		L[x...][d]
	elseif d ∈ -D:-1
		x = mod1.(addtuple(-1, -d, p), size(L))
		L[x...][-d]'
	else
		throw(ArgumentError("Link direction must be in range [1, $D] or [-$D, -1]. Got d = $d"))
	end
end
getlink(L::LocalLattice{D}, d::Int, p::Vararg{Int, D}) where D = getlink(L, d, Tuple(p))
getlink(L::LocalLattice{D}, d::Int, p::CartesianIndex{D}) where D = getlink(L, d, Tuple(p))

evenmask(A) = [iseven(sum(Tuple(i))) for i in CartesianIndices(A)]

end # module StaticQCD