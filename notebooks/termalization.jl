### A Pluto.jl notebook ###
# v0.19.3

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
using DistributedQCD, Plots, DataFrames, DataFramesMeta, PlutoUI

# ╔═╡ aa1a4f63-8df5-4ea3-a0dc-3ab59a0908f5
df = collect_results(datadir("polyloop"))

# ╔═╡ bb4b3a80-9310-44b1-ad13-2562fc6502bc
md"""
Observable to plot:
$(@bind col Select([
	:polyloop => "exp. val. of Polyakov loops", 
	:mod_polyloop => "exp. val. of modulus of Polyakov loops",
	:avg_plaq => "average plaquette",
]))
"""

# ╔═╡ 79f08bdc-94f1-4265-a6d2-0729a537553c
v, β = getproperty.(df.df, col), df.β

# ╔═╡ 29feb57e-6e0b-43bd-8790-1ed06849df69
stephist(v, bins = 100, label = ["β = $(β[1])" "β = $(β[2])" "β = $(β[3])"])

# ╔═╡ a3b5c6ec-e692-49bf-8d15-e1133d10bea3
plot(v[3])

# ╔═╡ Cell order:
# ╠═a7a2033e-19da-45ca-aee7-3821cff5f597
# ╠═263b0ba3-d191-475d-8183-3babab74affe
# ╠═0aafddce-0c54-4df5-bf3e-87d684782c87
# ╠═aa1a4f63-8df5-4ea3-a0dc-3ab59a0908f5
# ╟─bb4b3a80-9310-44b1-ad13-2562fc6502bc
# ╟─29feb57e-6e0b-43bd-8790-1ed06849df69
# ╠═79f08bdc-94f1-4265-a6d2-0729a537553c
# ╠═a3b5c6ec-e692-49bf-8d15-e1133d10bea3
