using DocStringExtensions, DrWatson, Statistics, DataFrames, Measurements, Plots, JLD2, LegibleLambdas, LsqFit, GLM, CSV
import Distributions: quantile, Normal, TDist
import Base.MathConstants: γ

jld2dir(args...) = DrWatson.projectdir("jld2", args...)

model4 = LegibleLambdas.@λ (x, p) -> p[2] .+ p[3] .* (x .- p[1]) .^ 2 + p[4] .* (x .- p[1]) .^ 3 + p[5] .* (x .- p[1]) .^ 4
model3 = LegibleLambdas.@λ (x, p) -> p[2] .+ p[3] .* (x .- p[1]) .^ 2 + p[4] .* (x .- p[1]) .^ 3
model2 = LegibleLambdas.@λ (x, p) -> p[2] .+ p[3] .* (x .- p[1]) .^ 2

model = model4 # default model

betamodel = LegibleLambdas.@λ (nt, p) -> p[1] .* nt .+ p[2] .+ p[3] ./ nt
linearmodel = LegibleLambdas.@λ (nt, p) -> p[1] .+ p[2] .* nt

# Bessel function
@. K₀(t) = sqrt(π / (2t)) * ℯ^(-t)

@. longdistance(R, p) = p[2] * (K₀(R / p[1]) + K₀((100 - R) / p[1]))

# Caristo's paper uses t = R / ξ
@. shortdistance(R, p) = p[2]/(R^(1/4)) * (
	1 + 
	R/(2p[1]) * log((ℯ^γ * R) / (8p[1])) +
	(R / p[1])^2 / 16 + 
	(R / p[1])^3 * log((ℯ^γ * R) / (8p[1])) / 32
)

convertkeys(d::Dict{String}) = Dict(Symbol.(keys(d)) .=> values(d))
convertkeys(d::Dict{Symbol}) = Dict(String.(keys(d)) .=> values(d))

mval(x::Measurements.Measurement) = getproperty(x, :val)
merr(x::Measurements.Measurement) = getproperty(x, :err)

mval(x) = mval.(x)
merr(x) = merr.(x)


### chi_squared ###

@doc "$(TYPEDSIGNATURES)"
function chi_squared(y, f)
	acc = 0.0
	for (yi, yr, fi) in zip(mval(y), merr(y), f)
		acc += (yi - fi)^2 / yr^2
	end
	acc
end

@doc "$(TYPEDSIGNATURES)"
function chi_squared(model::Function, fit::LsqFit.LsqFitResult, x, y)
	m = model(x, LsqFit.coef(fit))
	chi_squared(y, m) / LsqFit.dof(fit)
end

@doc "$(TYPEDSIGNATURES)"
chi_squared(fit::LsqFit.LsqFitResult, x, y) = chi_squared(model, fit, x, y)

@doc "$(TYPEDSIGNATURES)"
chi_squared(model::Function, fit::LsqFit.LsqFitResult, df::DataFrames.DataFrame) = chi_squared(model, fit, df.beta, df.susc)

@doc "$(TYPEDSIGNATURES)"
chi_squared(fit::LsqFit.LsqFitResult, df::DataFrames.DataFrame) = chi_squared(model, fit, df)

@doc "$(TYPEDSIGNATURES)"
function chi_squared(model::Function, d::Dict, nt::Integer)
	key = Symbol("nt$nt")
	df = d[:peak_points][key]
	chi_squared(model, d[:fit][key], df)
end

@doc "$(TYPEDSIGNATURES)"
chi_squared(d::Dict, nt::Integer) =	chi_squared(model, d, nt)


### fitmodel ###

"""
$(TYPEDSIGNATURES)

`y` must be a vector of Measurements.
"""
function fitmodel(model::Function, x, y, params; weighted = true, label1 = "y1", label2 = "y2", kwargs...)
	yavg = mval(y)

	fit = if weighted
		@info "Fitting with weights"
		w = 1 ./ merr(y) .^ 2
		LsqFit.curve_fit(model, x, yavg, w, params)
	else
		@info "Fitting without weights"
		LsqFit.curve_fit(model, x, yavg, params)
	end

	f(x) = model(x, LsqFit.coef(fit))
	m = f(x)
	chi = chi_squared(y, m) / LsqFit.dof(fit)

	p = plot(x, y; label = label1, kwargs...)
	plot!(p, f, label = label2)

	(fit = fit, plot = p, chi = chi, f = f)
end

@doc "$(TYPEDSIGNATURES)"
fitmodel(x, y, params; kwargs...) = fitmodel(model, x, y, params; kwargs...)

@doc "$(TYPEDSIGNATURES)"
fitmodel(model::Function, df::DataFrames.DataFrame, params; kwargs...) = fitmodel(model, df.beta, df.susc, params; kwargs...)

@doc "$(TYPEDSIGNATURES)"
fitmodel(df::DataFrames.DataFrame, params; kwargs...) = fitmodel(model, df, params; kwargs...)

@doc "$(TYPEDSIGNATURES)"
function fitmodel(model::Function, d::Dict, nt::Integer, params; kwargs...) 
	key = Symbol("nt$nt")
	fitmodel(model, d[:peak_points][key], params; kwargs...)
end

@doc "$(TYPEDSIGNATURES)"
fitmodel(d::Dict, nt::Integer, params; kwargs...) = fitmodel(model, d, nt, params; kwargs...)


### plot ###

@doc "$(TYPEDSIGNATURES)"
function Plots.plot(model::Function, d::Dict, nt::Integer; kwargs...)
	key = Symbol("nt$nt")
	x = d[:peak_points][key].beta
	y = d[:peak_points][key].susc
	f(t) = model(t, LsqFit.coef(d[:fit][key]))

	p = plot(x, y; kwargs...)
	plot!(p, f)
	p
end

@doc "$(TYPEDSIGNATURES)"
Plots.plot(d::Dict, nt::Integer; kwargs...) = Plots.plot(model, d, nt; kwargs...)

@doc "$(TYPEDSIGNATURES)"
function get_temperatures(linearfit, beta, nts)
	c = coef(linearfit)
	nt(β) = (β - c[1]) / c[2]
	T_c = 1 / nt(beta)
	T = 1 ./ nts
	T ./ T_c
end

### autocorrelation ###

function Γ(y, t; mean = nothing)
	m = if isnothing(mean)
		Statistics.mean(y)
	else
		mean
	end
	L = length(y)
	sum([(y[i] - m) * (y[i+t] - m) for i in 1:L-t]) / (L - t)
end

function τ(y, cutoff)
	m = Statistics.mean(y)
	d = Γ(y, 0, mean = m)
	(1 + 2*sum([Γ(y, t, mean = m) / d for t in 1:cutoff])) / 2
end