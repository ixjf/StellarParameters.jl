### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 0ba701f0-9006-11eb-08c9-6183e18fa170
begin
	using Pkg
	Pkg.activate("..")
	
	using Plots
	gr()
end

# ╔═╡ 25a5f9d0-9006-11eb-3c3b-eb70551516f0
include("../src/StellarParameters.jl")

# ╔═╡ 91aa6940-9006-11eb-357e-c5e62eee94cd
full_line_list, ew_list_Sun = read_line_list()

# ╔═╡ a7403aa0-9006-11eb-2825-3d0ac0de4494
multiplet_list = group_lines_into_multiplets(full_line_list)

# ╔═╡ af2bd120-9006-11eb-2eb3-87c224af2b6a
sisma_archive_path = joinpath(@__DIR__, "../../data/sisma_archive/")

# ╔═╡ b5562fa0-9006-11eb-0e28-cdf834acd9e7
begin
	file = joinpath(sisma_archive_path, "HD169822_20120724_0000_nor.fits")
	
	λ, flux, opt_data = read_spectrum_sisma(file)
	
	vrad = opt_data["VRAD"]
	
	lines, ew_list = isolate_all_lines_found_in_spectrum(
		vrad,
		nothing,
		keys(full_line_list),
		(λ, flux))
	
	grid_cache = GridCache()
	
	best_combo_i, best_combo_j, _ = find_best_multiplet_combo_for_Texc_estimate(
		weighted_median,
		multiplet_list, 
		Dict([λ_c_approx => (λ_c_approx,) for (λ_c_approx,_) in lines]), 
		ew_list_Sun)

	_, _, Texc = find_best_multiplet_combo_for_Texc_estimate(
		weighted_median,
		[multiplet_list[best_combo_i], multiplet_list[best_combo_j]],
		lines,
		ew_list)
	
	Teff, logg, FeH, α_Fe = match_ew_against_grid(ew_list, Texc, opt_data["R"], opt_data["VSINI"], grid_cache)
end

# ╔═╡ Cell order:
# ╠═0ba701f0-9006-11eb-08c9-6183e18fa170
# ╠═25a5f9d0-9006-11eb-3c3b-eb70551516f0
# ╠═91aa6940-9006-11eb-357e-c5e62eee94cd
# ╠═a7403aa0-9006-11eb-2825-3d0ac0de4494
# ╠═af2bd120-9006-11eb-2eb3-87c224af2b6a
# ╠═b5562fa0-9006-11eb-0e28-cdf834acd9e7
