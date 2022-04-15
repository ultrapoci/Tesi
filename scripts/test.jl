using DrWatson
@quickactivate "Tesi"

includet(srcdir("CabibboMarinari.jl"))

using QCD, ProgressMeter, Plots, Distributions

##

function rand_a()
	ϕ = rand(Uniform(0.0, 2π))
	θ = acos(rand(Uniform(-1.0, 1.0)))

	a₁ = sin(θ) * cos(ϕ)
	a₂ = sin(θ) * sin(ϕ)
	a₃ = cos(θ)

	a₁, a₂, a₃
end

A, B, C = [], [], []
for _ in 1:10^6
	(a, b, c) = rand_a()
	push!(A, a)
	push!(B, b)
	push!(C, c)
end