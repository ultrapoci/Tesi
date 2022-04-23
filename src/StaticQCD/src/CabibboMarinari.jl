export subrepresentations, generate_a0, randomSU2, overrelaxation!, heatbath!, touch_overrelaxation, touch_heatbath, normalizelattice!, one_termalization!

const subrepresentations = [ 
	(
		to_su2 = x::StaticMatrix{4, 4, ComplexF64} -> SU2(x[1, 1], x[1, 3]) |> normalizeSU2det,
		from_su2 = x::SU2 -> Sp2([x.t₁ 0; 0 1], [x.t₂ 0; 0 0])
	),
	(
		to_su2 = x::StaticMatrix{4, 4, ComplexF64} -> SU2(x[2, 2], x[2, 4]) |> normalizeSU2det,
		from_su2 = x::SU2 -> Sp2([1 0; 0 x.t₁], [0 0; 0 x.t₂])
	),
	(
		to_su2 = x::StaticMatrix{4, 4, ComplexF64} -> SU2(x[1, 1] + x[2, 2], x[1, 4] + x[2, 3]) |> normalizeSU2det,
		from_su2 = x::SU2 -> Sp2([x.t₁ 0; 0 x.t₁], [0 x.t₂; x.t₂ 0])
	),
	(
		to_su2 = x::StaticMatrix{4, 4, ComplexF64} -> SU2(x[1, 1] + x[4, 4], x[1, 2] - x[4, 3]) |> normalizeSU2det,
		from_su2 = x::SU2 -> Sp2([x.t₁ x.t₂; -conj(x.t₂) conj(x.t₁)], [0 0; 0 0])
	)
]

"""
	generate_a0(k::Real, β::Real)
Generates a real number ``a₀`` according to the distribution P(a₀) = √(1 - a₀^2) exp(a₀ β k).		
"""
function generate_a0(k::Real, β::Real)
	if k < 0
		throw(ArgumentError("k must be positive, got $k."))
	elseif β < 0
		throw(ArgumentError("β must be positive, got $β."))
	end

	reject = true
	a₀ = 0.0
	while reject
		x = rand(Uniform(exp(-2*β*k), 1.0)) #? are range extremes a problem?
		a₀ = 1 + log(x) / (β*k)
		reject = 1 - √(1 - a₀^2) > rand(Uniform(0.0, 1.0))
	end	
	a₀
end

function randomSU2(k::Real, β::Real)
	a₀ = generate_a0(k, β)
	ϕ = rand(Uniform(0.0, 2π))
	θ = acos(rand(Uniform(-1.0, 1.0)))

	r = √(1 - a₀^2)

	a₁ = r * sin(θ) * cos(ϕ)
	a₂ = r * sin(θ) * sin(ϕ)
	a₃ = r * cos(θ)

	SU2(complex(a₀, a₃), complex(a₂, a₁))
end

function touch_overrelaxation(S::Sp2, R::StaticMatrix{4, 4, ComplexF64})
	for (to_su2, from_su2) in subrepresentations
		a, _ = to_su2(S * R)
		S = from_su2(a^-2) * S
	end
	S
end

function touch_heatbath(S::Sp2, R::StaticMatrix{4, 4, ComplexF64}, β::Real)
	for (to_su2, from_su2) in subrepresentations
		a, k = to_su2(S * R)
		a = randomSU2(k, β) * a^-1
		S = from_su2(a) * S
	end
	S
end

"""
	overrelaxation!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, nover::Int; log = false) where D
	overrelaxation!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, nover::Int; kwargs...) where D
	overrelaxation!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, β::Real; kwargs...) where D
Do `nover` iterations of overrelaxation to the lattice `L`. `mask` is a boolean `DArray` with the same dimensions \
and distributed in the same way of `L`. The same is true for `inds`, only it contains the CartesianIndices of `L`.

The named tuple returned by the function `newlattice` can be passed directly as the first argument.
"""
function overrelaxation!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, nover::Int; log = false) where D
	for i in 1:nover, parity in (:even, :odd), u in 1:D
		log && @info "Overrelaxation" i parity u
		with_workers() do _
			mask = parity == :even ? evenmask[:L] : .!evenmask[:L]

			# filters indices using the mask `mask`
			for x in CartesianIndices(L[:L])[mask]
				link = L[:L][x][u]
				R = sumstaples(L, u, inds[:L][x])
				S = touch_overrelaxation(link, R)
				L[:L][x][u] = S
			end
		end
	end
end
overrelaxation!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, nover::Int; kwargs...) where D = overrelaxation!(T..., nover; kwargs...)
overrelaxation!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, nover::Int; kwargs...) where D = overrelaxation!(T..., nover; kwargs...)

"""
	heatbath!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, β::Real; log = false) where D
	heatbath!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, β::Real; kwargs...) where D 
	heatbath!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, β::Real; kwargs...) where D
Same as `overrelaxation`, only it applies the heat-bath algorithm to the lattice `L`, and does only one iteration.

`β` is the coefficient in front of the action in the exponential of the Boltzmann distribution.
"""
function heatbath!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, β::Real; log = false) where D
	for parity in (:even, :odd), u in 1:D
		log && @info "Heat bath" parity u
		with_workers() do _
			mask = parity == :even ? evenmask[:L] : .!evenmask[:L]

			# filters indices using the mask `mask`
			for x in CartesianIndices(L[:L])[mask]
				link = L[:L][x][u]
				R = sumstaples(L, u, inds[:L][x])
				S = touch_heatbath(link, R, β)
				L[:L][x][u] = S
			end
		end
	end
end
heatbath!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, β::Real; kwargs...) where D = heatbath!(T..., β; kwargs...)
heatbath!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, β::Real; kwargs...) where D = heatbath!(T..., β; kwargs...)

"""
	normalizelattice!(L::Lattice{D}; log = false) where D
	normalizelattice!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}; kwargs...) where D
	normalizelattice!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}; kwargs...) where D
Normalize all links in the lattice `L` to make sure they belong to ``Sp(2)`` and ``SU(4)``. It is possible to directly pass a tuple or named tuple containing the lattice. 
"""
function normalizelattice!(L::Lattice{D}; log = false) where D
	log && @info "Normalizing..."
	with_workers() do _
		for u in 1:D, x in eachindex(L[:L])
			L[:L][x][u] = normalizeSp2(L[:L][x][u])
		end
	end
end
normalizelattice!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}; kwargs...) where D = normalizelattice!(T[1]; kwargs...)
normalizelattice!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}; kwargs...) where D = normalizelattice!(T.lattice; kwargs...)

function one_termalization!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, nover::Int, β::Real, do_normalization = false; log = false) where D
	overrelaxation!(L, evenmask, inds, nover, log = log)
	heatbath!(L, evenmask, inds, β, log = log)
	do_normalization && normalizelattice!(L)
end
one_termalization!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, nover::Int, β::Real, do_normalization = false; kwargs...) where D = one_termalization!(T..., nover, β, do_normalization; kwargs...)
one_termalization!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, nover::Int, β::Real, do_normalization = false; kwargs...) where D = one_termalization!(T..., nover, β, do_normalization; kwargs...)