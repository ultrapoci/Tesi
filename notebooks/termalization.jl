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

# ╔═╡ a7a2033e-19da-45ca-aee7-3821cff5f597
using DrWatson

# ╔═╡ 263b0ba3-d191-475d-8183-3babab74affe
# ╠═╡ show_logs = false
@quickactivate "Tesi"

# ╔═╡ 0aafddce-0c54-4df5-bf3e-87d684782c87
using DistributedQCD, Plots, DataFrames, DataFramesMeta, PlutoUI, Statistics, Measurements, LaTeXStrings

# ╔═╡ 79fad70f-d410-47b3-9586-37bc96049add
include(srcdir("Utilities.jl"))

# ╔═╡ d5095d25-a5ed-432f-b5dc-8e5f2c68f35f
md"""
# Plots
"""

# ╔═╡ 35092cc8-7995-4af9-adcf-44738b9b7fc6
md"""
Folder to read: $(@bind folder Select(
(x->x=>x).(readdir(datadir("rust"))) |> collect
))
"""

# ╔═╡ aa1a4f63-8df5-4ea3-a0dc-3ab59a0908f5
results = collect_results(datadir(folder));

# ╔═╡ 9a0ec951-a1e0-4247-8497-4872cca3fe9e
md"""
X axis label: $(@bind xaxis_label TextField(default="β = 8/g²"))
Latex: $(@bind xlabel_latex CheckBox())
"""

# ╔═╡ 906aa83b-511d-4153-b823-9e4d9235a96d
md"""
Plot folder: $(@bind plot_folder TextField(default = folder))
"""

# ╔═╡ 3a4ed019-6af7-47c3-8614-8478007cf99f
md"""
β min = $(@bind min_β TextField())

β max = $(@bind max_β TextField())
"""

# ╔═╡ a1be3f0d-7ccc-4a0c-a40a-0bbb181a05c5
md"""
Error bars: $(@bind use_error_bars CheckBox(default=true))
"""

# ╔═╡ 3eecae2c-2775-4cdc-8430-9253fccef3d0
md"""
$(@bind save_button Button("Save plot"))
"""

# ╔═╡ 47761aa2-7222-444b-934d-e0593b2cecbd
md"""
## Data manipulation
"""

# ╔═╡ 4ccc646e-0602-4d6b-a0bf-d412261f9598
md"""
### Sort by β and dims
"""

# ╔═╡ 79f08bdc-94f1-4265-a6d2-0729a537553c
begin 
	select!(results, [:β, :dims, :nterm, :data])
	sort!(results, [:β, :dims, order(:nterm, rev=true)])
end;

# ╔═╡ 745d3199-d361-4083-95ff-251d2e034fa0
results2 = @transform(
	results,
	@byrow :Cᵥ = measurement(
		var((:data).S), 
		var((:data).S) * sqrt(2 / (length((:data).S) - 1))
	) / prod(:dims)
);

# ╔═╡ 9b8096f2-1b1b-47ba-b669-648b7e5f7d51
md"""
### Calculate mean and std error from data
"""

# ╔═╡ f95043e2-e15d-4d8b-983e-27d6e58acdf9
begin
	d = DataFrame()
	for data in results2.data
		append!(
			d, 
			Dict(
				name => measurement(col) # col is a Vector
				for (name, col) in zip(names(data), eachcol(data))
			)
		)
	end
	df = hcat(results2, d)
end;

# ╔═╡ 3794fe45-e2bf-4d49-9221-307687317b0e
md"""
### Specific heat
"""

# ╔═╡ b0e6fee8-c42b-45bf-a621-330051833314
obsnames = filter(x->x∉[:β, :dims, :nterm, :data], Symbol.(names(df)))

# ╔═╡ 813e6066-1a8f-41f3-b6bf-4482d3d8d8d9
md"""
Observable to plot: $(@bind plot_col Select(obsnames))
"""

# ╔═╡ edd17153-8b5c-4037-88ce-3eae67765268
md"""
Plot's title: $(@bind plot_title TextField(30, default=String(plot_col)))
Latex: $(@bind title_latex CheckBox())
"""

# ╔═╡ cf926d06-09f7-4845-8b6c-706e1d38c9bb
md"""
Y axis label: $(@bind yaxis_label TextField(default=String(plot_col)))
Latex: $(@bind ylabel_latex CheckBox())
"""

# ╔═╡ 2edfa35d-ae20-4fce-8ba3-4262fcd6f517
md"""
File name: $(@bind file_name TextField(default = string(plot_col)*".png"))
"""

# ╔═╡ 0b1bc5a8-e298-41b4-8d2e-648f3c6911e0
md"""
## Plotting $plot_col
"""

# ╔═╡ 7cb4536b-e3bd-4724-a50f-b82a4c0151ef
gf = groupby(df, :dims);

# ╔═╡ 66142c02-0d2c-4d16-bf6b-deb0b17f16fa
md"""
# Imported packages
"""

# ╔═╡ ba1a70d4-5eae-49ca-9837-966155cf2444
markers = [:circle, :rect, :diamond, :hexagon, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star5, :star6, :star7, :star8, :vline, :hline, :+, :x, :cross]

# ╔═╡ 43ceb641-52e5-4e89-bf2c-032b70d5c58f
function build_plot()
	p = plot(
		title = title_latex ? L"%$plot_title" : plot_title, 
		xlabel = xlabel_latex ? L"%$xaxis_label" : xaxis_label, 
		ylabel = ylabel_latex ? L"%$yaxis_label" : yaxis_label, 
		xminorticks = 5,
		yminorticks = 5, 
		legend = :outertopright,
		dpi = 300
	)
	for (g, marker) in zip(gf, Iterators.cycle(markers))
		L = g.dims[begin][2]
		v = use_error_bars ? getproperty(g, plot_col) : Measurements.value.(getproperty(g, plot_col))
		
		inds = try
			findall(x -> parse(Float64, min_β) ≤ x ≤ parse(Float64, max_β), g.β)
		catch
			1:length(g.β)
		end
		
		plot!(p, g.β[inds], v[inds], 
			label = "L = $L", 
			marker = (marker, 3.5), 
			msc = :auto,
			markerstrokewidth = 0.5
		)
	end	
	p
end;

# ╔═╡ ecdb29af-b963-48a8-8110-d4ad1a646b0a
p = build_plot()

# ╔═╡ 656191c0-2ae1-4797-9ca5-31e8425b84fb
md"""
## Saving variables
"""

# ╔═╡ 6305b71c-1722-4125-b276-c1ca93951ffb
save_plot = Ref(false);

# ╔═╡ 0098cd8f-4ac4-4fab-bb49-9271ea388e65
first_trigger = Ref(false);

# ╔═╡ 49008af8-d913-41fc-b990-b98a197fe1c9
begin
	save_button
	save_plot[] = true
	save_trigger = true
end;

# ╔═╡ 1212a794-8e23-423f-9f04-50d9be906314
begin
	save_trigger
	if first_trigger[] && save_plot[] 
		wsave(plotsdir(plot_folder, file_name), build_plot())
		@info "$file_name saved to $(plotsdir(plot_folder))"
	end
	save_plot[] = false
	first_trigger[] = true
end;

# ╔═╡ Cell order:
# ╟─d5095d25-a5ed-432f-b5dc-8e5f2c68f35f
# ╟─35092cc8-7995-4af9-adcf-44738b9b7fc6
# ╠═aa1a4f63-8df5-4ea3-a0dc-3ab59a0908f5
# ╟─813e6066-1a8f-41f3-b6bf-4482d3d8d8d9
# ╟─edd17153-8b5c-4037-88ce-3eae67765268
# ╟─9a0ec951-a1e0-4247-8497-4872cca3fe9e
# ╟─cf926d06-09f7-4845-8b6c-706e1d38c9bb
# ╟─906aa83b-511d-4153-b823-9e4d9235a96d
# ╟─2edfa35d-ae20-4fce-8ba3-4262fcd6f517
# ╟─3a4ed019-6af7-47c3-8614-8478007cf99f
# ╟─a1be3f0d-7ccc-4a0c-a40a-0bbb181a05c5
# ╟─3eecae2c-2775-4cdc-8430-9253fccef3d0
# ╟─1212a794-8e23-423f-9f04-50d9be906314
# ╟─0b1bc5a8-e298-41b4-8d2e-648f3c6911e0
# ╟─ecdb29af-b963-48a8-8110-d4ad1a646b0a
# ╠═43ceb641-52e5-4e89-bf2c-032b70d5c58f
# ╟─47761aa2-7222-444b-934d-e0593b2cecbd
# ╟─4ccc646e-0602-4d6b-a0bf-d412261f9598
# ╠═79f08bdc-94f1-4265-a6d2-0729a537553c
# ╠═745d3199-d361-4083-95ff-251d2e034fa0
# ╟─9b8096f2-1b1b-47ba-b669-648b7e5f7d51
# ╠═f95043e2-e15d-4d8b-983e-27d6e58acdf9
# ╟─3794fe45-e2bf-4d49-9221-307687317b0e
# ╠═b0e6fee8-c42b-45bf-a621-330051833314
# ╠═7cb4536b-e3bd-4724-a50f-b82a4c0151ef
# ╟─66142c02-0d2c-4d16-bf6b-deb0b17f16fa
# ╠═a7a2033e-19da-45ca-aee7-3821cff5f597
# ╠═263b0ba3-d191-475d-8183-3babab74affe
# ╠═0aafddce-0c54-4df5-bf3e-87d684782c87
# ╠═79fad70f-d410-47b3-9586-37bc96049add
# ╠═ba1a70d4-5eae-49ca-9837-966155cf2444
# ╟─656191c0-2ae1-4797-9ca5-31e8425b84fb
# ╠═6305b71c-1722-4125-b276-c1ca93951ffb
# ╠═0098cd8f-4ac4-4fab-bb49-9271ea388e65
# ╠═49008af8-d913-41fc-b990-b98a197fe1c9
