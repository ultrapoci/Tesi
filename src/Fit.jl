using DocStringExtensions, DrWatson, Statistics, DataFrames, Measurements, Plots, JLD2, LegibleLambdas, LsqFit, GLM, CSV, LaTeXStrings
import Distributions: quantile, Normal, TDist
import Base.MathConstants: γ
import SpecialFunctions

includet(srcdir("Jackknife.jl"))

picsdir(args...) = joinpath(raw"E:\Università\2020-2021\Tesi\tesi_doc\pics", args...)
desktopdir(args...) = joinpath(raw"C:\Users\Niky\Desktop", args...)
jld2dir(args...) = DrWatson.projectdir("jld2", args...)

model4 = LegibleLambdas.@λ (x, p) -> p[2] .+ p[3] .* (x .- p[1]) .^ 2 + p[4] .* (x .- p[1]) .^ 3 + p[5] .* (x .- p[1]) .^ 4
model3 = LegibleLambdas.@λ (x, p) -> p[2] .+ p[3] .* (x .- p[1]) .^ 2 + p[4] .* (x .- p[1]) .^ 3
model2 = LegibleLambdas.@λ (x, p) -> p[2] .+ p[3] .* (x .- p[1]) .^ 2

model = model4 # default model

betamodel = LegibleLambdas.@λ (nt, p) -> p[2] .* nt .+ p[1] .+ p[3] ./ nt
linearmodel = LegibleLambdas.@λ (nt, p) -> p[1] .+ p[2] .* nt

# Bessel function
@. K₀(t) = SpecialFunctions.besselk(0, t) # sqrt(π / (2t)) * ℯ^(-t)
@. K2₀(t) = sqrt(π / (2t)) * exp(-t) * (
	1
	- 1 / (8t)
	+ 9 / (2 * (8t)^2)
	- 225 / (6 * (8t)^3) 
)

@. longdistance80(R, p) = p[2] * (K₀(R / p[1]) + K₀((80 - R) / p[1]))
@. longdistance100(R, p) = p[2] * (K₀(R / p[1]) + K₀((100 - R) / p[1]))

@. longdistance80_2(R, p) = p[2] * (K2₀(R / p[1]) + K2₀((80 - R) / p[1]))
@. longdistance100_2(R, p) = p[2] * (K2₀(R / p[1]) + K2₀((100 - R) / p[1]))

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
function fitmodel(model::Function, x, y, params; weighted = true, label1 = "data points", label2 = "fit", kwargs...)
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
function temperature_from_nt(linearfit, beta, nts)
	c = measurement.(coef(linearfit), stderror(linearfit))
	nt(β) = (β - c[1]) / c[2]
	T_c = 1 / nt(beta)
	T = 1 ./ nts
	T ./ T_c
end

@doc "$(TYPEDSIGNATURES)"
function temperature_from_T(linearfit, beta, T)
	c = measurement.(coef(linearfit), stderror(linearfit))
	nt(β) = (β - c[1]) / c[2]
	T_c = 1 / nt(beta)
	T ./ T_c
end

@doc "$(TYPEDSIGNATURES)"
function temperature_from_beta(linearfit, nt, betas)
	c = measurement.(coef(linearfit), stderror(linearfit))
	β_c = c[1] + c[2] * nt
	betas ./ β_c
end

@doc "$(TYPEDSIGNATURES)"
function beta_from_temperature(linearfit, nt, T)
	c = measurement.(coef(linearfit), stderror(linearfit))
	beta_c = c[1] + c[2] * nt
	T * beta_c
end

@doc "$(TYPEDSIGNATURES)"
function beta_T_conv(linearfit, nt; beta = nothing, T = nothing)
	if !isnothing(beta)
		temperature_from_beta(linearfit, nt, beta)
	elseif !isnothing(T)
		beta_from_temperature(linearfit, nt, T)
	else
		throw(ArgumentError("Both beta and T are nothing"))
	end
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

function readfilelist(folder; betarange = nothing, betas = nothing, regex = r"beta=(.*).csv")
	filelist = readdir(folder, join=true)
	is_valid = if isnothing(betarange)
		if isnothing(betas)
			_ -> true
		else
			b -> b in betas
		end
	else
		b -> first(betarange) <= b <= last(betarange)
	end

	files = map(filelist) do f
		beta = parse(Float64, match(regex, f).captures[1])
		(beta, f)
	end

	files = filter(files) do (beta, _)
		is_valid(beta)
	end

	sort!(files, by = x -> first(x))

	map(files) do (beta, f)
		df = CSV.read(f, DataFrame)[200:end, :]
		(beta = beta, v = df.polyloop_sum_over_v)
	end
end

function readproffiles(folder)
	files = readdir(folder, join = true)
	map(files) do file
		filename = basename(file)
		@info "Reading $filename"

		info = split(filename, '_')
		l = parse(Int64, strip(info[2], 's'))
		nt = parse(Int64, strip(info[3], 't'))
		beta = parse(Float64, strip(info[4], 'b'))
		v = open(file) do f
			parse.(Float64, eachline(f))
		end

		(l = l, nt = nt, beta = beta, v = v)
	end
end

function montecarlohistory(rawdata, nt::Int, l::Int, beta::Real; range = nothing, kwargs...)
	i = findfirst(x -> x.beta == beta, rawdata["nt$nt"]["L$l"])
	if isnothing(i)
		throw("beta = $beta not found")
	else
		b = first(rawdata["nt$nt"]["L$l"][i])
		v = last(rawdata["nt$nt"]["L$l"][i])
		if isnothing(range)
			range = 1:length(v)
		end
		plot(
			v[range];
			dpi = 300,
			xlabel = "n",
			ylabel = "ϕ",
			legend = nothing,
			minorgrid = true,
			kwargs...
		)
	end
end

function makestephist(rawdata, nt::Int, l::Int, betafilter::Function; bins = 400, usetitle = false, kwargs...)
	polyloops = filter(x -> betafilter(x.beta), rawdata["nt$nt"]["L$l"])

	s = stephist(;
		dpi = 300,
		xlabel = "ϕ",
		ylabel = "P(ϕ)",
		legendtitle = "beta",
		kwargs...
	)

	for p in polyloops
		stephist!(s, p.v, bins = bins, label = p.beta)
	end

	if usetitle 
		plot!(s; title = "Polyakov loop, Nt = $nt, Ns = $l")
	end
	
	s
end

makestephist(rawdata, nt::Int, l::Int, (bmin, bmax)::Tuple{Real, Real}; kwargs...) = makestephist(rawdata, nt, l, b -> bmin <= b <= bmax; kwargs...)
makestephist(rawdata, nt::Int, l::Int, betas; kwargs...) = makestephist(rawdata, nt, l, b -> b in betas; kwargs...)

function makestephist(data, rawdata, nt::Int, l::Int; kwargs...)
	kn = "nt$nt"
	kl = "L$l"

	beta_c = mval(data[kn][kl]["beta_c"])
	v = rawdata[kn][kl]
	(_, minindex) = findmin(map(x -> abs(x.beta - beta_c), v))
	(lowerindex, higherindex) = clamp.((minindex - 2, minindex + 2), 1, length(v))
	betas = map(x -> x.beta, v[lowerindex:higherindex])
	makestephist(rawdata, nt, l, betas)
end

function getsusc(data, betarange = nothing; binsize = 5000, skip = nothing)
	isinrange = if isnothing(betarange)
		_ -> true
	else
		x -> first(betarange) <= x <= last(betarange)
	end

	points = filter(d -> isinrange(d.beta), data)
	susceptibility(ϕ², modϕ) = ϕ² - modϕ^2

	s = map(points) do p
		thermalized = if isnothing(skip)
			p.v
		else
			p.v[skip:end] # thermalization
		end
		phi2 = mean.(binsamples(thermalized .^ 2, binsize)) # uncorrelated data
		phimod = mean.(binsamples(abs.(thermalized), binsize)) # uncorrelated data
		y = jackknife(susceptibility, phi2, phimod)
		(beta = p.beta, y = y)
	end

	DataFrames.DataFrame(
		:beta => first.(s),
		:susc_over_v => mval.(last.(s)),
		:susc_over_v_error => merr.(last.(s)),
		:y => last.(s)
	)
end

function getsusc(data, nt, l, betarange = nothing; binsize = 5000, skip = nothing)
	isinrange = if isnothing(betarange)
		_ -> true
	else
		x -> first(betarange) <= x <= last(betarange)
	end

	points = sort(
		filter(d -> d.nt == nt && d.l == l && isinrange(d.beta), data), 
		by = d -> d.beta
	)

	susceptibility(ϕ², modϕ) = ϕ² - modϕ^2

	s = map(points) do p
		thermalized = if isnothing(skip)
			p.v
		else
			p.v[skip:end] # thermalization
		end
		phi2 = mean.(binsamples(thermalized .^ 2, binsize)) # uncorrelated data
		phimod = mean.(binsamples(abs.(thermalized), binsize)) # uncorrelated data
		y = jackknife(susceptibility, phi2, phimod)
		(beta = p.beta, y = y)
	end

	DataFrames.DataFrame(
		:beta => first.(s),
		:susc_over_v => mval.(last.(s)),
		:susc_over_v_error => merr.(last.(s)),
		:y => last.(s)
	)
end

function getsusc(rawdata::Dict, nt, l, betarange = nothing; binsize = 5000, skip = nothing)
	kn = "nt$nt"
	kl = "L$l"

	points = if isnothing(betarange)
		rawdata[kn][kl]
	else
		filter(
			p -> first(betarange) <= p.beta <= last(betarange), 
			rawdata[kn][kl]
		)
	end

	susceptibility(ϕ², modϕ) = ϕ² - modϕ^2

	s = map(points) do p
		thermalized = if isnothing(skip)
			p.v
		else
			p.v[skip:end] # thermalization
		end
		phi2 = mean.(binsamples(thermalized .^ 2, binsize)) # uncorrelated data
		phimod = mean.(binsamples(abs.(thermalized), binsize)) # uncorrelated data
		y = jackknife(susceptibility, phi2, phimod)
		(beta = p.beta, y = y)
	end

	DataFrames.DataFrame(
		:beta => first.(s),
		:susc_over_v => mval.(last.(s)),
		:susc_over_v_error => merr.(last.(s)),
		:y => last.(s)
	)
end

function getbinder(rawdata::Dict, nt, l, betarange = nothing; binsize = 5000, skip = nothing)
	kn = "nt$nt"
	kl = "L$l"

	points = if isnothing(betarange)
		rawdata[kn][kl]
	else
		filter(
			p -> first(betarange) <= p.beta <= last(betarange), 
			rawdata[kn][kl]
		)
	end

	binder(ϕ⁴, ϕ²) = ϕ⁴ / ϕ²^2 - 3.0

	s = map(points) do p
		thermalized = if isnothing(skip)
			p.v
		else
			p.v[skip:end] # thermalization
		end
		phi4 = mean.(binsamples(thermalized .^ 4, binsize)) # uncorrelated data
		phi2 = mean.(binsamples(thermalized .^ 2, binsize)) # uncorrelated data
		y = jackknife(binder, phi4, phi2)
		(beta = p.beta, y = y)
	end

	DataFrames.DataFrame(
		:beta => first.(s),
		:binder => mval.(last.(s)),
		:binder_error => merr.(last.(s)),
		:y => last.(s)
	)
end

function suscplot(data::Dict, nt, L; usetitle = false, kwargs...)
	kn = "nt$nt"
	kl = "L$L"
	df = data[kn][kl]["points"]

	f = if data[kn][kl]["model"] == "quartic"
		x -> model4(x, coef(data[kn][kl]["fit"]))
	elseif data[kn][kl]["model"] == "cubic"
		x -> model3(x, coef(data[kn][kl]["fit"]))
	else
		throw(ArgumentError("Found model = $(data[kn][kl]["model"])"))
	end

	chi = round(data[kn][kl]["r"].chi, digits = 4)
	beta_c = data[kn][kl]["beta_c"]
	m = data[kn][kl]["model"]
	binsize = data[kn][kl]["binsize"]
	p = scatter(
		df.beta,
		df.y;
		legend = :topleft,
		dpi = 300,
		xminorticks = 5,
		minorgrid = 5,
		titlefontsize = 10,
		xlabel = "beta",
		ylabel = "χ / Nₛ²",
		label = "data points",
		kwargs...
	)
	plot!(p, f, label = "fit ($m model)")

	if usetitle
		plot!(p, title = "Susceptibility fit: Nt=$nt, L=$L, binsize=$binsize\nβc=$beta_c, χ²=$chi")
	end

	p
end

#= function fssplot(data, kind::Symbol, nt, L = 60:20:100; kwargs...)
	kn = "nt$nt"

	(kk, exp_c, ylabel) = if kind == :modphi
		("modphi", 1/8, L"\langle|\phi|\rangle L^{\beta/\nu}")
	elseif kind == :phi2
		# -7/4 is gamma, 2 is to multiply by the volume
		("phi2", -7/4 + 2, L"L^2 \langle\phi^2\rangle L^{-\gamma/\nu}")
	else
		throw(ArgumentError("Curve's kind must be either :phi2 or :modphi, got :$kind"))
	end

	p = scatter(;
		dpi = 300,
		xlabel = L"xL^{1/\nu}",
		ylabel = ylabel,
		xminorticks = 5,
		minorgrid = true,
		legend = :topleft,
		title = "Nt = $nt",
		kwargs...
	)

	for l in L
		kl = "L$l"
		if haskey(data[kn], kl)
			beta_c = mval(data[kn][kl]["beta_c"])
			points = map(data[kn][kl][kk]) do (beta, y)
				x = beta / beta_c - 1
				(x * l, y * l^exp_c)
			end
			scatter!(p, 
				first.(points), 
				last.(points), 
				markershape = :cross,
				label = "L = $l"
			)
		else
			@info "Nt = $nt, L = $l not found, skipping"
		end
	end

	p
end =#

function fssplot(data, kind::Symbol, nt, L = 40:20:100; adj = nothing, usetitle = false, kwargs...)
	kn = "nt$nt"

	if adj == :auto
		return fssplot(data, kind, nt, L; adj = data[kn]["adj"][String(kind)], kwargs...)
	end

	(kk, exp_c, ylabel) = if kind == :modphi
		("modphi", 1/8, L"\langle|\phi|\rangle N_s^{\beta/\nu}")
	elseif kind == :phi2
		# -7/4 is gamma, 2 is to multiply by the volume
		("phi2", -7/4 + 2, L"L^2 \langle\phi^2\rangle N_s^{-\gamma/\nu}")
	else
		throw(ArgumentError("Curve's kind must be either :phi2 or :modphi, got :$kind"))
	end

	pl = scatter(;
		dpi = 300,
		xlabel = L"x N_s^{1/\nu}",
		ylabel = ylabel,
		xminorticks = 5,
		minorgrid = true,
		legend = :topleft,
		legendfontsize = 7,
		kwargs...
	)

	if usetitle
		plot!(pl, 
			title = if isnothing(adj)
				"Nt = $nt"
			elseif adj == :fit
				"Nt = $nt (beta_c adjusted fitting L=100)"
			else
				"Nt = $nt (beta_c adjusted)"
			end
		)
	end

	@. m(x, p) = p[3] / (1 + p[2] * exp(-p[1]*x))
	@. m_inv(y, p) = -log((p[3] / y - 1) / p[2]) / p[1]

	v = data[kn]["L100"][kk]
	beta_c100 = mval(data[kn]["L100"]["beta_c"])
	x = map(d -> 100 * (d.beta / beta_c100 - 1), v)
	y = map(d -> mval(last(d)) * 100^exp_c, v)
	fit = curve_fit(m, x, y, [1.0, 1.0, 1.0])
	f_inv(x) = m_inv(x, coef(fit))

	new_betas = []
	t_names = []
	t_points = []

	for (i, l) in enumerate(L)
		if l != 100
			kl = "L$l"
			if haskey(data[kn], kl)
				push!(t_names, Symbol(kl))
				if isnothing(adj)
					beta_c = mval(data[kn][kl]["beta_c"])
					points = map(data[kn][kl][kk]) do (beta, p)
						(x = (beta / beta_c - 1) * l, y = p * l^exp_c)
					end
					push!(new_betas, (L = l, beta = data[kn][kl]["beta_c"]))
					push!(t_points, points)
					scatter!(pl, 
						first.(points), 
						last.(points), 
						markershape = :cross,
						label = "Nₛ = $l"
					)
				elseif adj == :fit
					@info "Using fit to find new beta_c"

					beta_c = mval(data[kn][kl]["beta_c"])
					beta_err = merr(data[kn][kl]["beta_c"])
					v = map(data[kn][kl][kk]) do (beta, p)
						(beta = beta, x = l * (beta / beta_c - 1), y = p * l^exp_c)
					end

					(beta, _, ymin) = argmin(d -> abs(d.x), v)
					x2 = f_inv(mval(ymin))
					new_beta_c = beta * (x2 / l + 1)^(-1)

					push!(new_betas, (L = l, beta = measurement(new_beta_c, beta_err)))

					points = map(data[kn][kl][kk]) do (beta, p)
						(x = l * (beta / new_beta_c - 1), y = p * l^exp_c)
					end		

					push!(t_points, points)

					scatter!(pl, 
						first.(points), 
						last.(points), 
						markershape = :cross,
						label = "Nₛ = $l"
					)
				else
					beta_c = data[kn][kl]["beta_c"]
					new_beta_c = beta_c + adj[i]
					push!(new_betas, (L = l, beta = new_beta_c))

					points = map(data[kn][kl][kk]) do (beta, p)
						(x = l * (beta / mval(new_beta_c) - 1), y = p * l^exp_c)
					end

					push!(t_points, points)

					scatter!(pl, 
						first.(points), 
						last.(points), 
						markershape = :cross,
						label = "Nₛ = $l"
					)
				end
			else
				@info "Nt = $nt, L = $l not found, skipping"
			end
		end
	end

	if 100 in L
		push!(t_names, :L100)
		points = map(data[kn]["L100"][kk]) do (beta, p)
			(x = 100 * (beta / beta_c100 - 1), y = p * 100^exp_c)
		end
		push!(t_points, points)
		scatter!(pl, 
			first.(points), 
			last.(points), 
			markershape = :cross,
			label = "Nₛ = 100"
		)
		push!(new_betas, (L = 100, beta = data[kn]["L100"]["beta_c"]))
	end

	t = NamedTuple{tuple(t_names...)}(t_points)
	(plot = pl, betas = new_betas, adj = adj, t...)
end

function beta_vs_l(data, nt, L = 40:20:100; usetitle = false, kwargs...)
	kn = "nt$nt"
	b = map(L) do l 
		kl = "L$l"
		if haskey(data[kn], kl)
			(l, data[kn][kl]["beta_c"])
		else
			@info "Couldn't find nt=$nt L=$l, skipping"
			nothing
		end
	end
	filter!(x -> !isnothing(x), b)
	p = plot(
		first.(b),
		last.(b);
		dpi = 300,
		xlabel = L"N_s",
		ylabel = L"\beta_c",
		legend = nothing,
		kwargs...
	)

	if usetitle
		plot!(p; title = "Nt = $nt")
	end

	p
end

function getribbon(fit::LsqFit.LsqFitResult, xs, der; alpha = 0.05, dist::Symbol = :t)
	cov = LsqFit.estimate_covar(fit)
	D = if dist == :t
		TDist(dof(fit))
	elseif dist == :normal
		Normal()
	else
		throw(ArgumentError("'dist' must be either :t or :normal, got :$dist")) 
	end

	σ = map(xs) do x
		j = map(d -> d(x), der)
		sqrt(j' * cov * j)
	end

	z = quantile(D, 1 - alpha / 2)
	σ .* z
end

function betafit(data::Dict, model::Symbol; alpha = 0.05, dist::Symbol = :t, params = nothing, usetitle = false, kwargs...)
	(model, der, default_params, label, title) = if model == :linear
		(linearmodel, [_ -> 1, x -> x], [1.0, 1.0], "Linear fit", L"\beta_c(N_t) = a + bN_t")
	elseif model == :nonlinear
		(betamodel, [_ -> 1, x -> x, x -> 1 / x], [1.0, 1.0, 1.0], "Non linear fit", L"\beta_c(N_t) = a + bN_t + c/N_t")
	else
		throw(ArgumentError("'model' must be either :linear or :nonlinear, got :$model"))
	end

	if isnothing(params)
		params = default_params
	end

	nts = vcat(5:8...)
	betas = map(nt -> data["nt$nt"]["L100"]["beta_c"], nts)
	df = DataFrame(:nt => nts, :beta => betas)
	r = fitmodel(model, df.nt, df.beta, params) 
	
	c = measurement.(coef(r.fit), stderror(r.fit))
	c_names = ["a", "b", "c"]

	annotations = [((0.8, y), (rpad("$n = $k", 18), 10)) for (k, n, y) in zip(c, c_names, [0.5, 0.4, 0.3])]

	pl = plot(nts, r.f;
		ribbon = getribbon(r.fit, df.nt, der; alpha = alpha, dist = dist),
		minorgrid = true,
		xlabel = L"N_t",
		ylabel = L"\beta_c",
		dpi = 300,
		legend = :topleft,
		label = label,
		annotation = annotations,
		kwargs...
	)
	scatter!(pl, df.nt, df.beta, markershape = :cross, label = "Data points")

	if usetitle
		plot!(p, title = title * L",\ L=100,\ \chi^2 = %$(round(r.chi, digits = 4))")
	end

	(plot = pl, chi = r.chi, f = r.f, fit = r.fit, points = df)
end

function corrplot(data, nt, L; usetitle = false, kwargs...)
	dict = data["nt$nt"]["L$L"]["corr"]
	ks = keys(dict)

	
end