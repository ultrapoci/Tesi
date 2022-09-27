using DrWatson, Statistics, DataFrames, Measurements, Plots, JLD2, LegibleLambdas, LsqFit

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

function chi_squared(y, f)
	acc = 0.0
	for (yi, yr, fi) in zip(val(y), err(y), f)
		acc += (yi - fi)^2 / yr^2
	end
	acc
end

function chi_squared(model, fit, x, y)
	m = model(x, LsqFit.coef(fit))
	chi_squared(y, m) / LsqFit.dof(fit)
end

function fitmodel(model, df, params, weights = nothing; kwargs...)
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