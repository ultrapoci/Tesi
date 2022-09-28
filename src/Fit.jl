using DocStringExtensions, DrWatson, Statistics, DataFrames, Measurements, Plots, JLD2, LegibleLambdas, LsqFit

#= model = LegibleLambdas.@λ(
	(x, p) -> p[2] .+ p[3] .* (x .- p[1]) .^ 2 .+ p[4] .* (x .- p[1]) .^ 3 .+ p[5] .* (x .- p[1]) .^ 4
) =#

@. model(x, p) = p[2] + p[3] * (x - p[1])^2 + p[4] * (x - p[1])^3 + p[5] * (x - p[1])^4

convertkeys(d::Dict{String}) = Dict(Symbol.(keys(d)) .=> values(d))
convertkeys(d::Dict{Symbol}) = Dict(String.(keys(d)) .=> values(d))

val(x::Measurements.Measurement) = getproperty(x, :val)
err(x::Measurements.Measurement) = getproperty(x, :err)

val(x) = val.(x)
err(x) = err.(x)


### chi_squared ###

@doc "$(TYPEDSIGNATURES)"
function chi_squared(y, f)
	acc = 0.0
	for (yi, yr, fi) in zip(val(y), err(y), f)
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
function chi_squared(model::Function, d::Dict, nt::Integer)
	key = Symbol("nt$nt")
	df = d[:peak_points][key]
	chi_squared(model, d[:fit][key], df.beta, df.susc)
end

@doc "$(TYPEDSIGNATURES)"
chi_squared(d::Dict, nt::Integer) =	chi_squared(model, d, nt)


### fitmodel ###

@doc "$(TYPEDSIGNATURES)"
function fitmodel(model::Function, df::DataFrames.DataFrame, params, weights = nothing; kwargs...)
	xdata = df.beta
	ydata = df.susc_avg
	y = df.susc

	fit = if isnothing(weights) 
		@info "Fitting without weights"
		LsqFit.curve_fit(model, xdata, ydata, params)
	else
		@info "Fitting using weights"
		LsqFit.curve_fit(model, xdata, ydata, weights, params)
	end

	f(x) = model(x, LsqFit.coef(fit))
	m = f(xdata)
	chi = chi_squared(y, m) / LsqFit.dof(fit)

	p = plot(xdata, y; kwargs...)
	plot!(p, f)

	(fit = fit, plot = p, chi = chi, f = f)
end

@doc "$(TYPEDSIGNATURES)"
fitmodel(df::DataFrames.DataFrame, params, weights = nothing; kwargs...) = fitmodel(model, df, params, weights; kwargs...)

@doc "$(TYPEDSIGNATURES)"
function fitmodel(model::Function, d::Dict, nt::Integer, params, weights = nothing; kwargs...) 
	key = Symbol("nt$nt")
	fitmodel(model, d[:peak_points][key], params, weights; kwargs...)
end

@doc "$(TYPEDSIGNATURES)"
fitmodel(d::Dict, nt::Integer, params, weights = nothing; kwargs...) = fitmodel(model, d, nt, params, weights; kwargs...)


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