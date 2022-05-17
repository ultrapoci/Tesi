### A Pluto.jl notebook ###
# v0.19.4

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

# ╔═╡ b173768c-773f-4130-8029-64e2d2dd743d
md"""
## Read dataframe
"""

# ╔═╡ 35092cc8-7995-4af9-adcf-44738b9b7fc6
md"""
Folder to read: $(@bind folder Select(
(x->x=>x).(readdir(datadir())) |> collect
))
"""

# ╔═╡ aa1a4f63-8df5-4ea3-a0dc-3ab59a0908f5
results = collect_results(datadir(folder));

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

# ╔═╡ 9b8096f2-1b1b-47ba-b669-648b7e5f7d51
md"""
### Calculate mean and std error from data
"""

# ╔═╡ f95043e2-e15d-4d8b-983e-27d6e58acdf9
begin
	d = DataFrame()
	for data in results.data
		append!(
			d, 
			Dict(
				name => measurement(col) # col is a Vector
				for (name, col) in zip(names(data), eachcol(data))
			)
		)
	end
	df = hcat(results, d)
end;

# ╔═╡ 221b8424-5215-4e8f-a7ca-33c3b4fbe427
md"""
### Calculate susceptibility from <φ²> - <|φ|>²
"""

# ╔═╡ 50505008-1375-4d69-bea2-26f70fc703be

begin
	temp = DataFrame()
	for data in df.data
		V = 40*40
		append!(
			temp,
			Dict(:susc =>
				measurement(data.φ² .- 
				data.φ .^ 2)
			)
		)
	end
	data = hcat(df, temp);
end

#=
begin
	temp = DataFrame();
	for data in df.data
		append!(
			temp,
			Dict(:susc =>
				measurement(data.χᵥ .- data.polyloop .^ 2)
			)
		)
	end;
	data = hcat(df, temp);
end;
=#


#=data = @transform(df, :susc = 
	#measurement.(getproperty.(:data, ^(:φ²))) .- 
	#measurement.(getproperty.(:data, ^(:mod_φ))) .^ 2
	getproperty.(:data, ^(:φ²)) .- 
	(x->x.^2).(getproperty.(:data, ^(:mod_φ)))
);=#



#=data = @transform(df, :susc = 
	measurement.(
		getproperty.(:data, ^(:χᵥ)) .- 
		(x->x.^2).(getproperty.(:data, ^(:polyloop)))
	)
);=#


# ╔═╡ b0e6fee8-c42b-45bf-a621-330051833314
obsnames = filter(x->x∉[:β, :dims, :nterm, :data], Symbol.(names(data)))

# ╔═╡ 813e6066-1a8f-41f3-b6bf-4482d3d8d8d9
md"""
Observable to plot: $(@bind plot_col Select(obsnames))
"""

# ╔═╡ edd17153-8b5c-4037-88ce-3eae67765268
md"""
Plot's title: $(@bind plot_title TextField(default=String(plot_col)))

X axis label: $(@bind xaxis_label TextField(default="β = 8/g²"))

Y axis label: $(@bind yaxis_label TextField(default=String(plot_col)))

File name: $(@bind file_name TextField(default = "plot.png"))
"""

# ╔═╡ 0b1bc5a8-e298-41b4-8d2e-648f3c6911e0
md"""
## Plotting $plot_col
"""

# ╔═╡ 7cb4536b-e3bd-4724-a50f-b82a4c0151ef
gf = groupby(data, :dims);

# ╔═╡ 66142c02-0d2c-4d16-bf6b-deb0b17f16fa
md"""
# Imported packages
"""

# ╔═╡ ba1a70d4-5eae-49ca-9837-966155cf2444
markers = [:circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]

# ╔═╡ ecdb29af-b963-48a8-8110-d4ad1a646b0a
begin
	p = plot(
		title = plot_title, 
		xlabel = xaxis_label, 
		ylabel = yaxis_label, 
		xminorticks = 5,
		yminorticks = 5, 
		legend = :topleft,
		dpi = 300
	)
	for (g, marker) in zip(gf, Iterators.cycle(markers))
		L = g.dims[begin][2]
		v = getproperty(g, plot_col)
		plot!(p, g.β, v, 
			label = "L = $L", 
			marker = (marker, 3.5), 
			msc = :auto,
			markerstrokewidth = 0.5
		)
	end
	wsave(plotsdir("susc", file_name), p)
	@info "plot saved"
	p
end

# ╔═╡ Cell order:
# ╟─d5095d25-a5ed-432f-b5dc-8e5f2c68f35f
# ╟─813e6066-1a8f-41f3-b6bf-4482d3d8d8d9
# ╟─edd17153-8b5c-4037-88ce-3eae67765268
# ╟─0b1bc5a8-e298-41b4-8d2e-648f3c6911e0
# ╠═ecdb29af-b963-48a8-8110-d4ad1a646b0a
# ╟─b173768c-773f-4130-8029-64e2d2dd743d
# ╟─35092cc8-7995-4af9-adcf-44738b9b7fc6
# ╠═aa1a4f63-8df5-4ea3-a0dc-3ab59a0908f5
# ╟─47761aa2-7222-444b-934d-e0593b2cecbd
# ╟─4ccc646e-0602-4d6b-a0bf-d412261f9598
# ╠═79f08bdc-94f1-4265-a6d2-0729a537553c
# ╟─9b8096f2-1b1b-47ba-b669-648b7e5f7d51
# ╠═f95043e2-e15d-4d8b-983e-27d6e58acdf9
# ╟─221b8424-5215-4e8f-a7ca-33c3b4fbe427
# ╠═50505008-1375-4d69-bea2-26f70fc703be
# ╠═b0e6fee8-c42b-45bf-a621-330051833314
# ╠═7cb4536b-e3bd-4724-a50f-b82a4c0151ef
# ╟─66142c02-0d2c-4d16-bf6b-deb0b17f16fa
# ╠═a7a2033e-19da-45ca-aee7-3821cff5f597
# ╠═263b0ba3-d191-475d-8183-3babab74affe
# ╠═0aafddce-0c54-4df5-bf3e-87d684782c87
# ╠═79fad70f-d410-47b3-9586-37bc96049add
# ╠═ba1a70d4-5eae-49ca-9837-966155cf2444
