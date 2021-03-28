### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 40bbc780-8f17-11eb-343e-e1c079a00e28
begin
	using Pkg
	Pkg.activate("..")
	
	using Plots
	plotly()
	
	using Measurements: value, uncertainty
end

# ╔═╡ 652777f0-8f16-11eb-305e-33e87b6c651b
include("../src/StellarParameters.jl")

# ╔═╡ 91b36e00-8f16-11eb-2685-4fa4e7e7e892
full_line_list, ew_list_Sun = read_line_list()

# ╔═╡ 92aa5bb0-8f17-11eb-150d-db473910cdf6
multiplet_list = group_lines_into_multiplets(full_line_list)

# ╔═╡ c8d0f7d0-8f17-11eb-3bb0-83b0103721f1
sisma_archive_path = joinpath(@__DIR__, "../../data/sisma_archive/")

# ╔═╡ 8ac42a20-8f1c-11eb-355b-e1479c7b0dc5
begin
	d = Dict()
	
	for filename in filter(x -> endswith(x, "fits"), readdir(sisma_archive_path))
		λ, flux, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))
		
		vrad = opt_data["VRAD"]
		lines, ew_list = isolate_all_lines_found_in_spectrum(
			vrad,
			nothing,
			keys(full_line_list),
			(λ, flux))
		
		push!(d, filename => (λ, flux, opt_data, lines, ew_list))
	end
		
	d
end

# ╔═╡ feb4495e-8f17-11eb-3052-55640f8e9bfe
begin
	p = scatter([], [])

	Teffs = []
	weighted_multiplet_data = []
	unweighted_multiplet_data = []
	
	for (filename, (λ, flux, opt_data, lines, ew_list)) in d
		best_combo_i_weighted, best_combo_j_weighted, _ = find_best_multiplet_combo_for_Texc_estimate(
			weighted_median, 
			multiplet_list, 
			lines, 
			ew_list_Sun)

		_, _, Texc_weighted = find_best_multiplet_combo_for_Texc_estimate(
			weighted_median, 
			[multiplet_list[best_combo_i_weighted], multiplet_list[best_combo_j_weighted]],
			lines,
			ew_list)
		
		push!(weighted_multiplet_data, (best_combo_i_weighted, best_combo_j_weighted, Texc_weighted))
		
		best_combo_i_unweighted, best_combo_j_unweighted, _ = find_best_multiplet_combo_for_Texc_estimate(
			unweighted_median, 
			multiplet_list, 
			lines, 
			ew_list_Sun)

		_, _, Texc_unweighted = find_best_multiplet_combo_for_Texc_estimate(
			unweighted_median, 
			[multiplet_list[best_combo_i_unweighted], multiplet_list[best_combo_j_unweighted]],
			lines,
			ew_list)
		
		push!(unweighted_multiplet_data, (best_combo_i_unweighted, best_combo_j_unweighted, Texc_unweighted))

		Teff = opt_data["TEFF"]
		
		push!(Teffs, Teff)
	end

	scatter!(Teffs, map(x -> x[3], weighted_multiplet_data); label="Weighted median", markercolor=:royalblue2, markersize=2.5, markerstrokecolor=:royalblue2)
		scatter!(Teffs, map(x -> x[3], unweighted_multiplet_data); label="Unweighted median", markercolor=:grey22, markersize=2.5, markerstrokecolor=:grey22)
	
	plot!([minimum(Teffs), maximum(Teffs)], [minimum(Teffs), maximum(Teffs)]; label="Reference")
	xlabel!("Teff / K")
	ylabel!("Texc / K")
	
	p
end

# ╔═╡ Cell order:
# ╠═40bbc780-8f17-11eb-343e-e1c079a00e28
# ╠═652777f0-8f16-11eb-305e-33e87b6c651b
# ╠═91b36e00-8f16-11eb-2685-4fa4e7e7e892
# ╠═92aa5bb0-8f17-11eb-150d-db473910cdf6
# ╠═c8d0f7d0-8f17-11eb-3bb0-83b0103721f1
# ╠═8ac42a20-8f1c-11eb-355b-e1479c7b0dc5
# ╠═feb4495e-8f17-11eb-3052-55640f8e9bfe
