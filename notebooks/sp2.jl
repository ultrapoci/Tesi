### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 50226000-b84f-11ec-3345-abbbeab7adc1
begin
	using DrWatson
	DrWatson.@quickactivate "Tesi"
end

# ╔═╡ 0c972a75-b174-4f97-bb8e-1a6b36a0bdfd
using QCD

# ╔═╡ bc6cbc75-41e3-4104-8c07-ae25579dd499
l = Lattice(4, 4)

# ╔═╡ 7baf7be4-5363-44c6-b248-bed027cb8289
link = getlink(l, 1, 1, 1)

# ╔═╡ 7479f415-e2d3-4791-a63e-f61229afb9d8
getstaple(l, link, 2)

# ╔═╡ Cell order:
# ╠═50226000-b84f-11ec-3345-abbbeab7adc1
# ╠═0c972a75-b174-4f97-bb8e-1a6b36a0bdfd
# ╠═bc6cbc75-41e3-4104-8c07-ae25579dd499
# ╠═7baf7be4-5363-44c6-b248-bed027cb8289
# ╠═7479f415-e2d3-4791-a63e-f61229afb9d8
