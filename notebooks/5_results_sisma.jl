### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 61fd5c00-9f7b-11eb-19cf-6b9c7ed93f39
begin
	using Pkg
	Pkg.activate("..")
	
	using Plots
	gr()
	
	using Measurements: value, uncertainty, ±
	nothing
end

# ╔═╡ f735c3db-1694-451e-a2ec-a7a99a996d60
include("../src/StellarParameters.jl");

# ╔═╡ de5b857d-059a-408e-99fd-6952043d7cfd
full_line_list, ew_list_Sun = read_line_list()

# ╔═╡ eabb07d2-f9ea-40bd-8ebf-5fb7a9879c97
multiplet_list = group_lines_into_multiplets(full_line_list)

# ╔═╡ af32accb-7961-411f-8679-f1e8c4f45eef
sisma_archive_path = joinpath(@__DIR__, "../../data/sisma_archive/")

# ╔═╡ d82c9225-659a-476f-b92e-6a6b16c5dcc2
begin
	results = []
	
	grid_cache = GridCache()
	
	for filename in filter(x -> endswith(x, "fits"), readdir(sisma_archive_path))
		λ, flux, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))
		 
		lines, ew_list = isolate_all_lines_found_in_spectrum(
			opt_data["VRAD"],
			keys(full_line_list),
			(λ, flux),
			5, # already determined in 2_lines.jl
			5, # already determined in 2_lines.jl
			12.5) # already determined in 2_lines.jl
		
		_, _, Texc = find_best_multiplet_combo_for_Texc_estimate(
			weighted_median, 
			[multiplet_list[1], multiplet_list[4]],
			lines,
			ew_list,
			2)
		
		params = match_ew_against_grid(ew_list, Texc, grid_cache)
		
		isnothing(params) && continue

		Teff, logg, FeH, α_Fe = params

		push!(results, filename => (λ, flux, opt_data, Teff, logg, FeH, α_Fe))
	end
end

# ╔═╡ 0d4cfe5f-9b2c-4a6a-907c-299a49cf76e3
begin
	p_Teff = scatter([], [])
	p_logg = scatter([], [])
	p_FeH = scatter([], [])
	
	min_Teff = +Inf
	max_Teff = -Inf
	
	min_logg = +Inf
	max_logg = -Inf
	
	min_FeH = +Inf
	max_FeH = -Inf
	
	for (_, (_, _, opt_data, Teff, logg, FeH, _)) in results		
		Teff_ref = opt_data["TEFF"]
		logg_ref = opt_data["LOGG"]
		FeH_ref = opt_data["FEH"]
		
		global min_Teff = min(Teff_ref, min_Teff)
		global max_Teff = max(Teff_ref, max_Teff)
		
		global min_logg = min(logg_ref, min_logg)
		global max_logg = max(logg_ref, max_logg)
		
		global min_FeH = min(FeH_ref, min_FeH)
		global max_FeH = max(FeH_ref, max_FeH)
	end
	
	xlabel!(p_Teff, "Teff ref / K")
	ylabel!(p_Teff, "Teff obtained / K")
	plot!(p_Teff, [min_Teff, max_Teff], [min_Teff, max_Teff]; label="Reference", legend=nothing, ribbon=400) # SISMA uncertainty (highest)
	plot!(p_Teff, [min_Teff, max_Teff], [min_Teff, max_Teff]; label="Reference", legend=nothing, ribbon=240) # SISMA uncertainty

	xlabel!(p_logg, "logg ref / dex")
	ylabel!(p_logg, "logg obtained / dex")
	plot!(p_logg, [min_logg, max_logg], [min_logg, max_logg]; label="Reference", legend=nothing, ribbon=0.50) # SISMA uncertainty (highest)
	plot!(p_logg, [min_logg, max_logg], [min_logg, max_logg]; label="Reference", legend=nothing, ribbon=0.32) # SISMA uncertainty
	
	xlabel!(p_FeH, "[Fe/H] ref / dex")
	ylabel!(p_FeH, "[Fe/H] obtained / dex")
	plot!(p_FeH, [min_FeH, max_FeH], [min_FeH, max_FeH]; label="Reference", legend=nothing, ribbon=0.40) # SISMA uncertainty (highest)
	plot!(p_FeH, [min_FeH, max_FeH], [min_FeH, max_FeH]; label="Reference", legend=nothing, ribbon=0.33) # SISMA uncertainty
	
	for (_, (_, _, opt_data, Teff, logg, FeH, _)) in results
		Teff_ref = opt_data["TEFF"]
		logg_ref = opt_data["LOGG"]
		FeH_ref = opt_data["FEH"]
		
		scatter!(p_Teff, [Teff_ref], [(value(Teff) ± (uncertainty(Teff) + 250))]; label="Teff obtained")
		scatter!(p_logg, [logg_ref], [(value(logg) ± (uncertainty(logg) + 0.5))]; label="logg obtained")
		scatter!(p_FeH, [FeH_ref], [(value(FeH) ± (uncertainty(FeH) + 0.25))]; label="[Fe/H] obtained")
	end
end

# ╔═╡ 39f6e0cf-7078-4918-9ae1-a60fb901aa02
p_Teff

# ╔═╡ 0579e0d4-6b2b-49ac-8a6b-6086e115f58b
p_logg

# ╔═╡ af905b6d-c32b-45a0-a0b4-87c3fa4acdc7
p_FeH

# ╔═╡ f2f822e4-e384-478b-9f04-5cb4a51c3426
md"Teff obtained are mostly within the uncertainty range. There is an odd point out, but it is to be expected of this star (studies don't agree on its parameters). For [Fe/H], this is true as well. However, for logg, there seems to be a horizontal trend, rather than one in the direction of the reference line, implying that perhaps the algorithm is not determining logg correctly.

I should note that the reference parameters have an uncertainty associated to them as well. From the SISMA article, the mean uncertainties are ± 240 K for Teff, ± 0.32 dex for logg, and ± 0.33 dex for [Fe/H] for stars with Teff < 6000 K. These uncertainties can be almost twice as high for early-type stars: ± 400 K for Teff, ± 0.50 dex for logg, ± 0.40 dex for [Fe/H]."

# ╔═╡ Cell order:
# ╠═61fd5c00-9f7b-11eb-19cf-6b9c7ed93f39
# ╠═f735c3db-1694-451e-a2ec-a7a99a996d60
# ╠═de5b857d-059a-408e-99fd-6952043d7cfd
# ╠═eabb07d2-f9ea-40bd-8ebf-5fb7a9879c97
# ╠═af32accb-7961-411f-8679-f1e8c4f45eef
# ╠═d82c9225-659a-476f-b92e-6a6b16c5dcc2
# ╠═0d4cfe5f-9b2c-4a6a-907c-299a49cf76e3
# ╠═39f6e0cf-7078-4918-9ae1-a60fb901aa02
# ╠═0579e0d4-6b2b-49ac-8a6b-6086e115f58b
# ╠═af905b6d-c32b-45a0-a0b4-87c3fa4acdc7
# ╟─f2f822e4-e384-478b-9f04-5cb4a51c3426
