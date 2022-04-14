using LinearAlgebra

const id = ComplexF64.(
	[1 0 0
	 0 1 0
	 0 0 1]
)

const λ = map(x -> ComplexF64.(x), [
	[0 1 0
	 1 0 0
	 0 0 0],

	[0 -im 0
	 im 0  0
	 0  0  0],

	[1  0  0
	 0 -1  0
	 0  0  0],

	[0 0 1
	 0 0 0
	 1 0 0],

	[0  0 -im
	 0  0  0
	 im 0  0],

	[0 0 0
	 0 0 1
	 0 1 0],

	[0  0  0
	 0  0 -im
	 0 im  0],

	√3^-1 * 
	[1 0 0
	 0 1 0
	 0 0 -2], 
])

##

function randSU3()
	a = rand(Float64, 9)
	s = sum(abs2.(a))
	a = a ./ s
	V = [[a[1]*id]; im .* a[2:9] .* λ]
	U = sum(V)
	U / (det(U)^3)
end

function buildSU3(v::Vector{<:Real})
	v = ComplexF64.(v)
	[v[1]+im*(v[2]+v[3]/√3) v[4]+im*v[5] v[6]+im*v[7]
	-v[4]+im*v[5] v[1]+im*(-v[2]+v[3]/√3) v[8]+im*v[9]
	-v[6]+im*v[7] -v[8]+im*v[9] v[1]-2im*v[3]/√3]
end