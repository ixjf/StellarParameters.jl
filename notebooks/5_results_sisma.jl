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
			20) # already determined in 2_lines.jl
		
		_, _, Texc = find_best_multiplet_combo_for_Texc_estimate(
			weighted_median, 
			[multiplet_list[1], multiplet_list[4]],
			lines,
			ew_list,
			2)
		
		params = match_ew_against_grid(ew_list, Texc, grid_cache)

		Teff, logg, FeH, α_Fe = params

		push!(results, filename => (λ, flux, opt_data, Teff, logg, FeH, α_Fe))
	end
end

# ╔═╡ 0d4cfe5f-9b2c-4a6a-907c-299a49cf76e3
begin
	p_Teff = scatter([], []; label="")
	p_logg = scatter([], []; label="")
	p_FeH = scatter([], []; label="")
	
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
	
	xlabel!(p_Teff, "Nº do espetro")
	ylabel!(p_Teff, "Teff obtido - ref / K")
	yticks!(p_Teff, -2500:250:1000)
	hline!(p_Teff, [0]; label="Intervalo incerteza máximo", ribbon=400) # SISMA uncertainty (highest)
	hline!(p_Teff, [0]; label="Intervalo incerteza", ribbon=240, legend=:bottomright) # SISMA uncertainty
	#plot!(p_Teff, [min_Teff, max_Teff], [min_Teff, max_Teff]; label="Reference", legend=nothing, ribbon=400) # SISMA uncertainty (highest)
	#plot!(p_Teff, [min_Teff, max_Teff], [min_Teff, max_Teff]; label="Reference", legend=nothing, ribbon=240) # SISMA uncertainty

	xlabel!(p_logg, "Nº do espetro")
	ylabel!(p_logg, "logg obtido - ref / dex")
	yticks!(p_logg, -2.5:0.5:3)
	hline!(p_logg, [0]; label="Intervalo incerteza máximo", ribbon=0.50) # SISMA uncertainty (highest)
	hline!(p_logg, [0]; label="Intervalo incerteza", ribbon=0.32, legend=:topright) # SISMA uncertainty
	#plot!(p_logg, [min_logg, max_logg], [min_logg, max_logg]; label="Reference", legend=nothing, ribbon=0.50) # SISMA uncertainty (highest)
	#plot!(p_logg, [min_logg, max_logg], [min_logg, max_logg]; label="Reference", legend=nothing, ribbon=0.32) # SISMA uncertainty
	
	xlabel!(p_FeH, "Nº do espetro")
	ylabel!(p_FeH, "[Fe/H] obtido - ref / dex")
	yticks!(p_FeH, -1.5:0.25:1.75)
	hline!(p_FeH, [0]; label="Intervalo incerteza máximo", ribbon=0.40) # SISMA uncertainty (highest)
	hline!(p_FeH, [0]; label="Intervalo incerteza", ribbon=0.33, legend=:topright) # SISMA uncertainty
	#plot!(p_FeH, [min_FeH, max_FeH], [min_FeH, max_FeH]; label="Reference", legend=nothing, ribbon=0.40) # SISMA uncertainty (highest)
	#plot!(p_FeH, [min_FeH, max_FeH], [min_FeH, max_FeH]; label="Reference", legend=nothing, ribbon=0.33) # SISMA uncertainty
	
	i = 0
	for (_, (_, _, opt_data, Teff, logg, FeH, _)) in results
		global i += 1
		
		Teff_ref = opt_data["TEFF"]
		logg_ref = opt_data["LOGG"]
		FeH_ref = opt_data["FEH"]
		
		scatter!(p_Teff, [i], [Teff - Teff_ref]; label="", markercolor=:grey)
		scatter!(p_logg, [i], [logg - logg_ref]; label="", markercolor=:grey)
		scatter!(p_FeH, [i], [FeH - FeH_ref]; label="", markercolor=:grey)
	end
end

# ╔═╡ 39f6e0cf-7078-4918-9ae1-a60fb901aa02
p_Teff

# ╔═╡ 60282e0d-4aec-43a5-a1c5-dac79461959e
savefig(p_Teff, "resultados_sisma_Teff.pdf")

# ╔═╡ 0579e0d4-6b2b-49ac-8a6b-6086e115f58b
p_logg

# ╔═╡ b80dbbc8-f378-434a-9501-cc4c2e562e69
savefig(p_logg, "resultados_sisma_logg.pdf")

# ╔═╡ af905b6d-c32b-45a0-a0b4-87c3fa4acdc7
p_FeH

# ╔═╡ 06b35a86-fd20-4106-9e73-086384b864a2
savefig(p_FeH, "resultados_sisma_FeH.pdf")

# ╔═╡ f2f822e4-e384-478b-9f04-5cb4a51c3426
md"Teff obtained are mostly within the uncertainty range. There is an odd point out, but it is to be expected of this star (studies don't agree on its parameters). For [Fe/H], this is true as well. However, for logg, there seems to be a horizontal trend, rather than one in the direction of the reference line, implying that perhaps the algorithm is not determining logg correctly.

I should note that the reference parameters have an uncertainty associated to them as well. From the SISMA article, the mean uncertainties are ± 240 K for Teff, ± 0.32 dex for logg, and ± 0.33 dex for [Fe/H] for stars with Teff < 6000 K. These uncertainties can be almost twice as high for early-type stars: ± 400 K for Teff, ± 0.50 dex for logg, ± 0.40 dex for [Fe/H]."

# ╔═╡ 99bd1539-87c4-4814-ab2b-5ee94d52f409
md"No graph for [α/Fe] since there are no reference values for this parameter for the SISMA spectra."

# ╔═╡ Cell order:
# ╠═61fd5c00-9f7b-11eb-19cf-6b9c7ed93f39
# ╠═f735c3db-1694-451e-a2ec-a7a99a996d60
# ╠═de5b857d-059a-408e-99fd-6952043d7cfd
# ╠═eabb07d2-f9ea-40bd-8ebf-5fb7a9879c97
# ╠═af32accb-7961-411f-8679-f1e8c4f45eef
# ╠═d82c9225-659a-476f-b92e-6a6b16c5dcc2
# ╠═0d4cfe5f-9b2c-4a6a-907c-299a49cf76e3
# ╠═39f6e0cf-7078-4918-9ae1-a60fb901aa02
# ╠═60282e0d-4aec-43a5-a1c5-dac79461959e
# ╠═0579e0d4-6b2b-49ac-8a6b-6086e115f58b
# ╠═b80dbbc8-f378-434a-9501-cc4c2e562e69
# ╠═af905b6d-c32b-45a0-a0b4-87c3fa4acdc7
# ╠═06b35a86-fd20-4106-9e73-086384b864a2
# ╟─f2f822e4-e384-478b-9f04-5cb4a51c3426
# ╟─99bd1539-87c4-4814-ab2b-5ee94d52f409
