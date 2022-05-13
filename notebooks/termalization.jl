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

# ╔═╡ d5095d25-a5ed-432f-b5dc-8e5f2c68f35f
md"""
# Plots
"""

# ╔═╡ 47761aa2-7222-444b-934d-e0593b2cecbd
md"""
## Data manipulation
"""

# ╔═╡ 4ccc646e-0602-4d6b-a0bf-d412261f9598
md"""
### Sort by β and dims
"""

# ╔═╡ ad47f215-ba50-45d6-b105-c584e30c85cc
md"""
### Remove duplicates
"""

# ╔═╡ e2d78392-435c-4001-8596-f42da7858e07
md"""
### Subtract the condensate from χ/L³
"""

# ╔═╡ dff10576-2feb-418c-903e-ef11e47f8d74
md"""
### Group data by dims
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
df = collect_results(datadir(folder));

# ╔═╡ 79f08bdc-94f1-4265-a6d2-0729a537553c
begin 
	select!(df, [:β, :dims, :nterm, :data])
	sort!(df, [:β, :dims, order(:nterm, rev=true)])
end;

# ╔═╡ 8d29521a-f371-4307-bcda-465900aa235d
@subset!(df, 
	.!((:β .== 6.46).&(:dims .∈ [[(2,8,8,8), (2,10,10,10), (2,12,12,12)]]).&(:nterm .== 300)));

# ╔═╡ dcf4af44-188d-428d-9d4f-6ae10923c429
data = 
@transform(df, 
	:χᵥ = (x->measurement(mean(x), std(x))).(
		getproperty.(:data, ^(:χᵥ)) .- 
		(x->x.^2).(getproperty.(:data, ^(:polyloop)))
	)
);

# ╔═╡ 7cb4536b-e3bd-4724-a50f-b82a4c0151ef
gf = groupby(data, :dims);

# ╔═╡ 66142c02-0d2c-4d16-bf6b-deb0b17f16fa
md"""
# Imported packages
"""

# ╔═╡ ba1a70d4-5eae-49ca-9837-966155cf2444
markers = [:circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]

# ╔═╡ a4622874-cec8-486a-a669-24f6db8db510
begin
	p = plot(
		title = L"\chi / L^3 = \langle\phi\rangle^2 -\langle\phi^2\rangle" * " for (3+1)d Sp(2)", 
		xlabel = "β = 8/g²", 
		ylabel = "χ/L³", 
		xminorticks = 5,
		yminorticks = 5, 
		legend = :bottomleft,
		dpi = 300
	)
	for (g, marker) in zip(gf, Iterators.cycle(markers))
		L = g.dims[begin][2]
		plot!(g.β, g.χᵥ, 
			label = "L = $L", 
			marker = (marker, 3.5), 
			msc = :auto,
			markerstrokewidth = 0.2
		)
	end
	wsave(plotsdir("susc", "fig8.png"), p)
	@info "plot saved"
	p
end

# ╔═╡ Cell order:
# ╟─d5095d25-a5ed-432f-b5dc-8e5f2c68f35f
# ╟─a4622874-cec8-486a-a669-24f6db8db510
# ╟─47761aa2-7222-444b-934d-e0593b2cecbd
# ╟─4ccc646e-0602-4d6b-a0bf-d412261f9598
# ╠═79f08bdc-94f1-4265-a6d2-0729a537553c
# ╟─ad47f215-ba50-45d6-b105-c584e30c85cc
# ╠═8d29521a-f371-4307-bcda-465900aa235d
# ╟─e2d78392-435c-4001-8596-f42da7858e07
# ╠═dcf4af44-188d-428d-9d4f-6ae10923c429
# ╟─dff10576-2feb-418c-903e-ef11e47f8d74
# ╠═7cb4536b-e3bd-4724-a50f-b82a4c0151ef
# ╟─b173768c-773f-4130-8029-64e2d2dd743d
# ╟─35092cc8-7995-4af9-adcf-44738b9b7fc6
# ╠═aa1a4f63-8df5-4ea3-a0dc-3ab59a0908f5
# ╟─66142c02-0d2c-4d16-bf6b-deb0b17f16fa
# ╠═a7a2033e-19da-45ca-aee7-3821cff5f597
# ╠═263b0ba3-d191-475d-8183-3babab74affe
# ╠═0aafddce-0c54-4df5-bf3e-87d684782c87
# ╠═ba1a70d4-5eae-49ca-9837-966155cf2444
