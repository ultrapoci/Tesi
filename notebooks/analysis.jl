### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c035fbc0-f62c-11ec-3347-b7f3f6c2064b
using DrWatson

# ╔═╡ 6027f926-06c0-4301-bea7-c006864b011f
@quickactivate "Tesi"

# ╔═╡ 04e43061-35a8-4e43-b415-1cae0bd9226a
begin
	using DistributedQCD, Plots, Statistics, Measurements, CSV, DataFrames, PlutoUI, LaTeXStrings
end

# ╔═╡ 87534882-cf45-424e-a16b-b1f2a206dc4c
include(srcdir("Utilities.jl"))

# ╔═╡ 90967b98-7cde-42a5-acf4-6a33add82be7
include(srcdir("Jackknife.jl"));

# ╔═╡ 9d70e646-4dfb-4fed-9b00-383fdf3923bb
md"""
# Plot
"""

# ╔═╡ 07a9828f-3de7-44e7-af62-4aa9ccac9cda
datafolder = datadir("rust", "nt=2")

# ╔═╡ 5d983982-21c0-4c83-b14b-a8e9fd6c9683
skip_offset = 400

# ╔═╡ 37178a96-1173-43c6-ba4a-9e5e8ee4ad32
binsize = 20

# ╔═╡ 8daccf89-a676-4772-ad1a-e9f74e7b9ccd
md"""
X axis label: $(@bind xaxis_label TextField(default="β = 8/g²"))
Latex: $(@bind xlabel_latex CheckBox())
"""

# ╔═╡ d9e369ce-bce3-4f64-bb2d-fe715769af1c
md"""
Plot folder: $(@bind plot_folder TextField())
"""

# ╔═╡ 53030e1a-f9cc-4c40-9545-fec41b099cd8
md"""
β min = $(@bind min_β TextField())

β max = $(@bind max_β TextField())
"""

# ╔═╡ bb7314f0-80c7-4e0d-89e9-2302824e40aa
md"""
Select L:

All: $(@bind all_ls CheckBox(default=true))
"""

# ╔═╡ f6d1d0c2-6ce6-4ff8-8ac7-3de03545ce16
md"""
Error bars: $(@bind use_error_bars CheckBox(default=true))
"""

# ╔═╡ f4afedf4-2f2d-465d-9218-4a227d3828d9
md"""
$(@bind save_button Button("Save plot"))
"""

# ╔═╡ 1fc2eae0-b2ce-4839-b2b7-658d2970352d
md"""
# Data manipulation
"""

# ╔═╡ 459314f0-efed-4231-8a68-fca7bd0da737
filelist = readdir(datafolder, join=true)

# ╔═╡ 8cdb6d2d-2248-472f-b883-64c33a940187
function readinfo(f)
	d = CSV.read(f, DataFrame, limit = 1)[1, :]
	Dict(Symbol.(names(d)) .=> values(d))
end

# ╔═╡ 74653267-7842-469a-a9d1-c5cd693cdec6
readdata(f) = CSV.read(f, DataFrame, header = 3)

# ╔═╡ 35a679fe-989b-4b27-950d-afd56840ae88
function add_measurements!(dict, offset = 0)
	data = last(dict[:data], nrow(dict[:data]) - offset)
	mapcols!(col -> mean.(binsamples(col, binsize)), data)
	
	V = dict[:v] # volume
	Vs = dict[:v_s] #spatial volume
	np = dict[:dims] * (dict[:dims] - 1) ÷ 2 # n° of plaquettes
	β = dict[:beta]
	
	avg_plaq(p) = p / (np * V) # average plaquette
	χ(ϕ², modϕ) = ϕ² - modϕ^2 / Vs # susceptibility
	χᵥ(ϕ², modϕ) = χ(ϕ², modϕ) / Vs # susceptibility per volume
	S(p) = -β * p / 4 # action
	S²(p) = S(p)^2 # action squared
	gᵣ(ϕ², ϕ⁴) = Vs * ϕ⁴ / ϕ²^2 - 3 # Binder cumulant
	Cᵥ(s², s) = (s² - s^2) / V # specific heat

	dict[:avg_plaq] = jackknife(avg_plaq, data.plaquettes)
	dict[:χ] = jackknife(χ, data.polyloops2, data.polyloops_mod)
	dict[:χᵥ] = jackknife(χᵥ, data.polyloops2, data.polyloops_mod)
	dict[:S] = jackknife(S, data.plaquettes)
	dict[:S²] = jackknife(S², data.plaquettes)
	dict[:gᵣ] = jackknife(gᵣ, data.polyloops2, data.polyloops4)
	dict[:Cᵥ] = jackknife(Cᵥ, S².(data.plaquettes), S.(data.plaquettes))
	#dict[:Cᵥ] = (mean(S².(data.plaquettes)) - mean(S.(data.plaquettes))) / V
	dict[:susc] = mean(χᵥ.(data.polyloops2, data.polyloops_mod))
	dict[:susc2] = χᵥ(mean(data.polyloops2), mean(data.polyloops_mod))
end

# ╔═╡ 37679ab2-f768-47ec-927d-d359d6fd036c
begin
	dicts = []
	for filename in filelist
		@info "Reading $filename"
		dict = readinfo(filename)	
		dict[:data] = readdata(filename)
		push!(dicts, dict)		
	end
end;

# ╔═╡ 31bc8726-638b-402e-8ff9-34920a452d31
begin
	add_measurements!.(dicts, skip_offset)
	tmpdf = DataFrame()
	for d in dicts
		append!(tmpdf, DataFrame(d))
	end
	df = sort(tmpdf, [:l, :beta])
end;

# ╔═╡ f4b0c054-7500-4fed-8fba-ab9b481d410d
dfs = groupby(df, :l);

# ╔═╡ 30ca6749-22b3-4134-88eb-48ccc2b79f0d
if !all_ls
	md"""$(@bind l_list MultiSelect(first.(keys(dfs)), default = first.(keys(dfs))))"""
end

# ╔═╡ 1d135b91-8121-43e3-86a9-e583a947ef32
obsnames = [:avg_plaq, :χ, :χᵥ, :S, :S², :gᵣ, :Cᵥ, :susc, :susc2]

# ╔═╡ e5fd0de8-8222-4e0d-a555-6cb078960e77
md"""
Observable to plot: $(@bind plot_col Select(obsnames))
"""

# ╔═╡ 2cab4f1c-9175-4696-b5f7-11a244acb67b
md"""
Plot's title: $(@bind plot_title TextField(30, default=String(plot_col)))
Latex: $(@bind title_latex CheckBox())
"""

# ╔═╡ ff191399-7cd3-408d-9487-9b96ef8f4b69
md"""
Y axis label: $(@bind yaxis_label TextField(default=String(plot_col)))
Latex: $(@bind ylabel_latex CheckBox())
"""

# ╔═╡ d2ed4a87-00d0-4270-a3ad-eb3331b7b317
md"""
File name: $(@bind file_name TextField(default = string(plot_col)*".png"))
"""

# ╔═╡ 4cdaf72c-57da-41c7-ac54-c0b6290d6830
markers = [:circle, :rect, :diamond, :hexagon, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star5, :star6, :star7, :star8, :vline, :hline, :+, :x, :cross]

# ╔═╡ f3ef5f6d-5646-4daa-95ad-6b84d14bf17e
begin
	p = plot(
		title = title_latex ? L"%$plot_title" : plot_title, 
		xlabel = xlabel_latex ? L"%$xaxis_label" : xaxis_label, 
		ylabel = ylabel_latex ? L"%$yaxis_label" : yaxis_label, 
		xminorticks = 5,
		yminorticks = 5, 
		legend = :outertopright,
		dpi = 300
	)

	Ls = try
		l_list
	catch
		first.(keys(dfs))
	end
	
	for (df, L, marker) in zip(dfs, first.(keys(dfs)), Iterators.cycle(markers))
		L ∉ Ls && continue
			
		v = use_error_bars ? getproperty(df, plot_col) : Measurements.value.(getproperty(df, plot_col))
		
		inds = try
			findall(x -> parse(Float64, min_β) ≤ x ≤ parse(Float64, max_β), df.beta)
		catch
			1:length(df.beta)
		end
		
		plot!(p, df.beta[inds], v[inds], 
			label = "L = $L", 
			marker = (marker, 3.5), 
			msc = :auto,
			markerstrokewidth = 0.5
		)
	end	
	p
end

# ╔═╡ 5085b1f1-581c-4f2a-a815-caa3d9f8f301
md"""
# Packages
"""

# ╔═╡ 6b265df7-cf51-4cb7-9846-c94471525308
md"""
# Save variables
"""

# ╔═╡ b573b0a5-a23a-47e0-9b6b-e69f49563e89
save_plot = Ref(false);

# ╔═╡ 19ffdd0b-1b22-4ee7-bc52-5042b9058280
first_trigger = Ref(false);

# ╔═╡ 0fd26021-cc9f-4383-bf32-43b11e564f0f
begin
	save_button
	save_plot[] = true
	save_trigger = true
end;

# ╔═╡ ea26e086-00d4-4996-b259-5e01caaa69e5
begin
	save_trigger
	if first_trigger[] && save_plot[] 
		wsave(plotsdir(plot_folder, file_name), p)
		@info "$file_name saved to $(plotsdir(plot_folder))"
	end
	save_plot[] = false
	first_trigger[] = true
end;

# ╔═╡ 2b6e4e88-db1e-46d4-b9b8-29a0b10cf615
md"""
# Test
"""

# ╔═╡ 4a430497-5767-4218-94fa-0cbdfac448d5
filename = filelist[1]

# ╔═╡ 09b55a16-789b-4bc0-be4d-712f933e7017
info = readinfo(filename)

# ╔═╡ 0ba2e0fb-1a69-44f0-967e-3d2ed4c423f4
test = readdata(filename)

# ╔═╡ 763db54e-2694-4e05-8574-9d6918fb4327
plaqs = test.plaquettes

# ╔═╡ a539170d-6cba-4e8c-9247-8bead91b42f5
mean(plaqs) / (3 * 8^3)

# ╔═╡ Cell order:
# ╟─9d70e646-4dfb-4fed-9b00-383fdf3923bb
# ╠═07a9828f-3de7-44e7-af62-4aa9ccac9cda
# ╠═5d983982-21c0-4c83-b14b-a8e9fd6c9683
# ╠═37178a96-1173-43c6-ba4a-9e5e8ee4ad32
# ╟─2cab4f1c-9175-4696-b5f7-11a244acb67b
# ╟─8daccf89-a676-4772-ad1a-e9f74e7b9ccd
# ╟─ff191399-7cd3-408d-9487-9b96ef8f4b69
# ╟─d9e369ce-bce3-4f64-bb2d-fe715769af1c
# ╟─d2ed4a87-00d0-4270-a3ad-eb3331b7b317
# ╟─53030e1a-f9cc-4c40-9545-fec41b099cd8
# ╟─bb7314f0-80c7-4e0d-89e9-2302824e40aa
# ╟─30ca6749-22b3-4134-88eb-48ccc2b79f0d
# ╟─f6d1d0c2-6ce6-4ff8-8ac7-3de03545ce16
# ╟─f4afedf4-2f2d-465d-9218-4a227d3828d9
# ╟─ea26e086-00d4-4996-b259-5e01caaa69e5
# ╟─e5fd0de8-8222-4e0d-a555-6cb078960e77
# ╠═f3ef5f6d-5646-4daa-95ad-6b84d14bf17e
# ╟─1fc2eae0-b2ce-4839-b2b7-658d2970352d
# ╟─459314f0-efed-4231-8a68-fca7bd0da737
# ╠═8cdb6d2d-2248-472f-b883-64c33a940187
# ╠═74653267-7842-469a-a9d1-c5cd693cdec6
# ╠═35a679fe-989b-4b27-950d-afd56840ae88
# ╠═37679ab2-f768-47ec-927d-d359d6fd036c
# ╠═31bc8726-638b-402e-8ff9-34920a452d31
# ╠═f4b0c054-7500-4fed-8fba-ab9b481d410d
# ╠═1d135b91-8121-43e3-86a9-e583a947ef32
# ╠═4cdaf72c-57da-41c7-ac54-c0b6290d6830
# ╟─5085b1f1-581c-4f2a-a815-caa3d9f8f301
# ╠═c035fbc0-f62c-11ec-3347-b7f3f6c2064b
# ╠═6027f926-06c0-4301-bea7-c006864b011f
# ╠═04e43061-35a8-4e43-b415-1cae0bd9226a
# ╠═87534882-cf45-424e-a16b-b1f2a206dc4c
# ╠═90967b98-7cde-42a5-acf4-6a33add82be7
# ╟─6b265df7-cf51-4cb7-9846-c94471525308
# ╠═b573b0a5-a23a-47e0-9b6b-e69f49563e89
# ╠═19ffdd0b-1b22-4ee7-bc52-5042b9058280
# ╠═0fd26021-cc9f-4383-bf32-43b11e564f0f
# ╟─2b6e4e88-db1e-46d4-b9b8-29a0b10cf615
# ╠═4a430497-5767-4218-94fa-0cbdfac448d5
# ╠═09b55a16-789b-4bc0-be4d-712f933e7017
# ╠═0ba2e0fb-1a69-44f0-967e-3d2ed4c423f4
# ╠═763db54e-2694-4e05-8574-9d6918fb4327
# ╠═a539170d-6cba-4e8c-9247-8bead91b42f5
