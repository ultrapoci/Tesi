const subrepresentations = [ 
	(
		to_su2 = x::StaticMatrix{4, 4, ComplexF64} -> SU2(x[1, 1], x[1, 4]) |> normalizeSU2det,
		from_su2 = x::SU2 -> Sp2([x.t₁ 0; 0 1], [0 x.t₂; 0 0])
	),
	(
		to_su2 = x::StaticMatrix{4, 4, ComplexF64} -> SU2(x[2, 2], x[2, 3]) |> normalizeSU2det,
		from_su2 = x::SU2 -> Sp2([1 0; 0 x.t₁], [0 0; x.t₂ 0])
	),
	(
		to_su2 = x::StaticMatrix{4, 4, ComplexF64} -> SU2(x[1, 1] + x[2, 2], x[1, 3] - x[2, 4]) |> normalizeSU2det,
		from_su2 = x::SU2 -> Sp2([x.t₁ 0; 0 x.t₁], [x.t₂ 0; 0 -x.t₂])
	),
	(
		to_su2 = x::StaticMatrix{4, 4, ComplexF64} -> SU2(x[1, 1] + x[3, 3], x[1, 2] + x[3, 4]) |> normalizeSU2det,
		from_su2 = x::SU2 -> Sp2([x.t₁ x.t₂; -conj(x.t₂) conj(x.t₁)], [0 0; 0 0])
	)
]

"""
	generate_a0(k::Real, β::Real)
Generates a real number ``a₀`` according to the distribution P(a₀) = √(1 - a₀^2) exp(a₀ β k).		
"""
function generate_a0(k::Real, β::Real)
	k < 0 && throw(ArgumentError("k must be positive, got $k."))
	β < 0 && throw(ArgumentError("β must be positive, got $β."))
	
	reject = true
	a₀ = 0.0
	while reject
		x = rand(Uniform(exp(-2*β*k), 1.0)) #? are range extremes a problem?
		a₀ = 1 + log(x) / (β*k)
		reject = (1 - a₀^2) < rand(Uniform(0.0, 1.0))^2
	end	
	a₀
end

function randomSU2(k::Real, β::Real)::SU2
	a₀ = generate_a0(k, β)
	ϕ = rand(Uniform(0.0, 2π))
	θ = acos(rand(Uniform(-1.0, 1.0)))

	r = √(1 - a₀^2)

	a₁ = r * sin(θ) * cos(ϕ)
	a₂ = r * sin(θ) * sin(ϕ)
	a₃ = r * cos(θ)

	SU2(complex(a₀, a₃), complex(a₂, a₁))
end

function touch_overrelaxation(S::Sp2, R::SMatrix{4, 4, ComplexF64})::Sp2
	U::Sp2 = S
	for (to_su2, from_su2) in subrepresentations
		a::SU2, _ = to_su2(U * R)::Tuple{SU2, Real}
		U = from_su2(a^-2) * U
	end
	U
end

function touch_heatbath(S::Sp2, R::SMatrix{4, 4, ComplexF64}, β::Real)::Sp2
	U::Sp2 = S
	for (to_su2, from_su2) in subrepresentations
		a::SU2, k::Real = to_su2(U * R)::Tuple{SU2, Real}
		a = randomSU2(k, β) * a^-1
		U = from_su2(a) * U
	end
	U
end


"""
	overrelaxation!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, nover::Int; log = false, iter = missing) where D
	overrelaxation!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, nover::Int; kwargs...) where D
	overrelaxation!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, β::Real; kwargs...) where D
Do `nover` iterations of overrelaxation to the lattice `L`. `mask` is a boolean `DArray` with the same dimensions \
and distributed in the same way of `L`. The same is true for `inds`, only it contains the CartesianIndices of `L`.

The named tuple returned by the function `newlattice` can be passed directly as the first argument.
"""
function overrelaxation!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, nover::Int; log = false, iter = missing) where D
	for i in 1:nover, parity in (:even, :odd), u in 1:D
		log && @info "Overrelaxation" iter cycle=i parity direction=u
		
		with_workers(procs = vec(procs(L))) do
			mask = parity == :even ? evenmask[:L] : .!evenmask[:L]

			# filters indices using the mask `mask`
			@maybe_threaded for x in CartesianIndices(L[:L])[mask] # x is the position of the link in the local array
				link::Sp2 = L[:L][x][u]
				R = sumstaples(L, u, inds[:L][x]) # inds[:L][x] is the position of the link in the global array 
				L[:L][x][u] = touch_overrelaxation(link, R)
			end
		end
	end
end
overrelaxation!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, nover::Int; kwargs...) where D = overrelaxation!(T..., nover; kwargs...)
overrelaxation!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, nover::Int; kwargs...) where D = overrelaxation!(T..., nover; kwargs...)

"""
	heatbath!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, β::Real; log = false, iter = missing) where D
	heatbath!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, β::Real; kwargs...) where D 
	heatbath!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, β::Real; kwargs...) where D
Same as `overrelaxation`, only it applies the heat-bath algorithm to the lattice `L`, and does only one iteration.

`β` is the coefficient in front of the action in the exponential of the Boltzmann distribution.
"""
function heatbath!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, β::Real; log = false, iter = missing) where D
	for parity in (:even, :odd), u in 1:D
		log && @info "Heat-bath" iter parity direction=u

		with_workers(procs = vec(procs(L))) do
			mask = parity == :even ? evenmask[:L] : .!evenmask[:L]

			# filters indices using the mask `mask`
			@maybe_threaded for x in CartesianIndices(L[:L])[mask] # x is the position of the link in the local array
				link::Sp2 = L[:L][x][u]
				R = sumstaples(L, u, inds[:L][x]) # inds[:L][x] is the position of the link in the global array 
				L[:L][x][u] = touch_heatbath(link, R, β)
			end
		end
	end
end
heatbath!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, β::Real; kwargs...) where D = heatbath!(T..., β; kwargs...)
heatbath!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, β::Real; kwargs...) where D = heatbath!(T..., β; kwargs...)

"""
	normalizelattice!(L::Lattice{D}; log = false, iter = missing) where D
	normalizelattice!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}; kwargs...) where D
	normalizelattice!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}; kwargs...) where D
Normalize all links in the lattice `L` to make sure they belong to ``\\mathrm{Sp(2)}`` and ``\\mathrm{SU(4)}``. It is possible to directly pass a tuple or named tuple containing the lattice. 
"""
function normalizelattice!(L::Lattice{D}; log = false, iter = missing) where D
	log && @info "Normalizing..." iter

	with_workers(procs = vec(procs(L))) do
		@maybe_threaded for u in 1:D
			for x in eachindex(L[:L])
				L[:L][x][u] = normalizeSp2(L[:L][x][u])
			end
		end
	end
end
normalizelattice!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}; kwargs...) where D = normalizelattice!(T[1]; kwargs...)
normalizelattice!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}; kwargs...) where D = normalizelattice!(T.lattice; kwargs...)

"""
	one_termalization!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, nover::Int, β::Real, do_normalization = false; log = false, iter = missing) where D
	one_termalization!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, nover::Int, β::Real, do_normalization = false; kwargs...) where D
	one_termalization!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, nover::Int, β::Real, do_normalization = false; kwargs...) where D
Do one iteration of termalization, meaning `nover` cycles of overrelaxation, one cycle of heatbath and, if `do_normalization` is true, normalize all lattice.
"""
function one_termalization!(L::Lattice{D}, evenmask::Mask{D}, inds::Indices{D}, nover::Int, β::Real, do_normalization = false; log = false, iter = missing) where D
	overrelaxation!(L, evenmask, inds, nover, log = log, iter = iter)
	heatbath!(L, evenmask, inds, β/2, log = log, iter = iter) # β/2 is to account for trace normalization
	do_normalization && normalizelattice!(L, log = log, iter = iter)
end
one_termalization!(T::Tuple{Lattice{D}, Mask{D}, Indices{D}}, nover::Int, β::Real, do_normalization = false; kwargs...) where D = one_termalization!(T..., nover, β, do_normalization; kwargs...)
one_termalization!(T::NamedTuple{(:lattice, :mask, :inds), Tuple{Lattice{D}, Mask{D}, Indices{D}}}, nover::Int, β::Real, do_normalization = false; kwargs...) where D = one_termalization!(T..., nover, β, do_normalization; kwargs...)