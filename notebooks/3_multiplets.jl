### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 40bbc780-8f17-11eb-343e-e1c079a00e28
begin
	using Pkg
	Pkg.activate("..")
	
	using Plots, StatsPlots
	gr()
	
	using Measurements: value, uncertainty, Measurement
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
			keys(full_line_list),
			(λ, flux),
			5, # already determined in 2_lines.jl
			5, # already determined in 2_lines.jl
			12.5) # already determined in 2_lines.jl
		
		push!(d, filename => (λ, flux, opt_data, lines, ew_list))
	end
end

# ╔═╡ feb4495e-8f17-11eb-3052-55640f8e9bfe
begin
	Teffs = []
	weighted_multiplet_data = []
	unweighted_multiplet_data = []
	
	for (filename, (λ, flux, opt_data, lines, ew_list)) in d
		best_combo_i_weighted, best_combo_j_weighted, _ = find_best_multiplet_combo_for_Texc_estimate(
			weighted_median, 
			multiplet_list, 
			Dict([λ_c_approx => (λ_c_approx,) for (λ_c_approx,_) in lines]), 
			ew_list_Sun,
			2)

		_, _, Texc_weighted = find_best_multiplet_combo_for_Texc_estimate(
			weighted_median, 
			[multiplet_list[best_combo_i_weighted], multiplet_list[best_combo_j_weighted]],
			lines,
			ew_list,
			2)
		
		push!(weighted_multiplet_data, (best_combo_i_weighted, best_combo_j_weighted, Texc_weighted))
		
		best_combo_i_unweighted, best_combo_j_unweighted, _ = find_best_multiplet_combo_for_Texc_estimate(
			unweighted_median, 
			multiplet_list, 
			Dict([λ_c_approx => (λ_c_approx,) for (λ_c_approx,_) in lines]), 
			ew_list_Sun,
			2)

		_, _, Texc_unweighted = find_best_multiplet_combo_for_Texc_estimate(
			unweighted_median, 
			[multiplet_list[best_combo_i_unweighted], multiplet_list[best_combo_j_unweighted]],
			lines,
			ew_list,
			2)
		
		push!(unweighted_multiplet_data, (best_combo_i_unweighted, best_combo_j_unweighted, Texc_unweighted))

		Teff = opt_data["TEFF"]
		
		push!(Teffs, Teff)
	end
end

# ╔═╡ 5cb65bae-8ff0-11eb-1d7f-79ea6ef8cf6b
begin
	p = scatter([], []; label="")
	
	scatter!(Teffs, map(x -> x[3], unweighted_multiplet_data); label="Unweighted median", markercolor=:grey22, markersize=2.5, markerstrokecolor=:grey22)
	scatter!(Teffs, map(x -> x[3], weighted_multiplet_data); label="Weighted median", markercolor=:blue1, markersize=2.5, markerstrokecolor=:blue1, markershape=:xcross, ylim=(3000, 8000))
	
	plot!([minimum(Teffs), maximum(Teffs)], [minimum(Teffs), maximum(Teffs)]; label="Reference")
	xlabel!("Teff / K")
	ylabel!("Texc / K")
	
	p
end

# ╔═╡ 26f7da10-9d3f-11eb-2bb4-1b32c3ec0d46
md"Both methods give fairly similar temperature estimates. The temperatures seem to work best around ~ 5750 K (temperature of the Sun), and then perhaps get worse and worse as the temperature of the star increases (however, notice how it produces an accurate result for the star with Texc ~ 4200 K, which is in accordance to other papers I've found). However, we can see that the error bars for the unweighted median are much larger than for the weighted median. Therefore, the weighted median should be chosen."

# ╔═╡ 87b410a0-8ff0-11eb-24b9-b7f48d7569b5
begin
	histo = Dict()
	for i=1:length(MULTIPLETS),j=i+1:length(MULTIPLETS)
		push!(histo, (i,j) => 1)
	end
	for (i,j,_) in unweighted_multiplet_data
		histo[(i,j)] += 1
	end
	#MULTIPLETS, ["($(x[1]), $(x[2]))" for x in (collect(keys(histo)))], reshape(collect(values(histo)), 1, length(histo)), repeat([""], inner=length(histo)) 
	groupedbar(["($(x[1]), $(x[2]))" for x in (collect(keys(histo)))], Int64.(reshape(collect(values(histo)), 1, length(histo))), group = repeat([""], inner=length(histo)), ylabel="No. times chosen", xlabel="Multiplet combo")
end

# ╔═╡ b257efa0-9d3f-11eb-1b79-53bfb7068707
md"One multiplet combo seems to be preferred. However, there are other combos which are also selected many times. I wonder if the difference in results is very large, or if they're chosen for a very small difference in the error during calibration."

# ╔═╡ 790eaae0-9eac-11eb-11d2-fdb7db273a7d
begin
	Teffs_1_4_combo = Float64[]
	multiplet_data_1_4_combo = Measurement{Float64}[]
	
	for (filename, (λ, flux, opt_data, lines, ew_list)) in d
		_, _, Texc = find_best_multiplet_combo_for_Texc_estimate(
			weighted_median, 
			[multiplet_list[1], multiplet_list[4]],
			lines,
			ew_list,
			2)
		
		push!(multiplet_data_1_4_combo, Texc)
		
		Teff = opt_data["TEFF"]
		
		push!(Teffs_1_4_combo, Teff)
	end
end

# ╔═╡ 815f1dae-9eac-11eb-0c31-fd90b64c8b44
begin
	p_1_4_combo = scatter([], []; label="")
	
	scatter!(Teffs, map(x -> x[3], weighted_multiplet_data); label="Texc (per-spectrum combo)", markercolor=:royalblue2, markersize=2.5, markerstrokecolor=:royalblue2)
	scatter!(p_1_4_combo, Teffs_1_4_combo, multiplet_data_1_4_combo; label="Texc (combo 1+4)", markercolor=:grey22, markersize=2.5, markerstrokecolor=:grey22)
	
	plot!([minimum(Teffs_1_4_combo), maximum(Teffs_1_4_combo)], [minimum(Teffs_1_4_combo), maximum(Teffs_1_4_combo)]; label="Reference")
	xlabel!("Teff / K")
	ylabel!("Texc / K")
	
	p_1_4_combo
end

# ╔═╡ 5d836650-9eaf-11eb-2ba3-c5d8ddd7d636
md"It turns out that while using the best combo per-spectrum gives the most accurate values for Texc, the uncertainties are actually significantly higher in many cases than when using combo 1+4. Since all uncertainties fall within range of the reference values, it seems best to always use combo 1+4."

# ╔═╡ 2e51291c-206e-45b3-b09d-ac84bd1c6ca9
md"Weighted median is used.

σ = 2

Previously determined values for filtering."

# ╔═╡ Cell order:
# ╠═40bbc780-8f17-11eb-343e-e1c079a00e28
# ╠═652777f0-8f16-11eb-305e-33e87b6c651b
# ╠═91b36e00-8f16-11eb-2685-4fa4e7e7e892
# ╠═92aa5bb0-8f17-11eb-150d-db473910cdf6
# ╠═c8d0f7d0-8f17-11eb-3bb0-83b0103721f1
# ╠═8ac42a20-8f1c-11eb-355b-e1479c7b0dc5
# ╠═feb4495e-8f17-11eb-3052-55640f8e9bfe
# ╠═5cb65bae-8ff0-11eb-1d7f-79ea6ef8cf6b
# ╟─26f7da10-9d3f-11eb-2bb4-1b32c3ec0d46
# ╠═87b410a0-8ff0-11eb-24b9-b7f48d7569b5
# ╟─b257efa0-9d3f-11eb-1b79-53bfb7068707
# ╠═790eaae0-9eac-11eb-11d2-fdb7db273a7d
# ╠═815f1dae-9eac-11eb-0c31-fd90b64c8b44
# ╟─5d836650-9eaf-11eb-2ba3-c5d8ddd7d636
# ╟─2e51291c-206e-45b3-b09d-ac84bd1c6ca9
