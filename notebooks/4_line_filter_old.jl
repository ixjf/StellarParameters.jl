### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ b2becdce-9eaa-11eb-3e6e-0fb4203c73ec
begin
	using Pkg
	Pkg.activate("..")
	
	using Plots, StatsPlots
	gr()
	
	using Measurements: value, uncertainty
end

# ╔═╡ c3c8ae70-9eaa-11eb-2d7d-310fcee9b3ed
include("../src/StellarParameters.jl")

# ╔═╡ c7670ae0-9eaa-11eb-3c60-6b9f4b0ce567
full_line_list, ew_list_Sun = read_line_list()

# ╔═╡ ca61fb60-9eaa-11eb-2222-dd88d37de4b7
multiplet_list = group_lines_into_multiplets(full_line_list)

# ╔═╡ ce2d8250-9eaa-11eb-067e-e19e370f8b05
sisma_archive_path = joinpath(@__DIR__, "../../data/sisma_archive/")

# ╔═╡ d1580e50-9eaa-11eb-1364-55020bd1fe06
begin
	all_spectra_ew_filter_comp = Dict()
	
	for filename in filter(x -> endswith(x, "fits"), readdir(sisma_archive_path))[1:20]
		ew_filter_comp = Dict()

		λ, flux, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))

		vrad = opt_data["VRAD"]

		grid_cache = GridCache()

		for N=0:5
			lines, ew_list = isolate_all_lines_found_in_spectrum(
				vrad,
				keys(full_line_list),
				(λ, flux),
				N,
				1,
				25)

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

			push!(ew_filter_comp, N => (Teff, logg, FeH, α_Fe))
		end
		
		push!(all_spectra_ew_filter_comp, filename => ew_filter_comp)
	end

	all_spectra_ew_filter_comp
end

# ╔═╡ ddee8f40-9eaa-11eb-23c9-ff7ebe92c868
begin
	p_Teff_ew_filter = scatter([], [], label="", legend=nothing)#:outerright)
	p_logg_ew_filter = scatter([], [], label="", legend=nothing)
	p_FeH_ew_filter = scatter([], [], label="", legend=nothing)
	#p_α_Fe_ew_filter = scatter([], [], label="", legend=nothing)
	#k = collect(keys(all_spectra_ew_filter_comp))
	#ew_filter_comp = all_spectra_ew_filter_comp[k[14]]
	
	for (filename, ew_filter_comp) in all_spectra_ew_filter_comp
		_, _, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))

		Teff_ref = opt_data["TEFF"]
		logg_ref = opt_data["LOGG"]
		FeH_ref = opt_data["FEH"]
			
		ew_filter_comp = sort!(collect(ew_filter_comp); by=x -> x[1])
		
		Ns = [N for (N, _) in ew_filter_comp]
		Teffs = [Teff for (_, (Teff, _, _, _)) in ew_filter_comp]
		loggs = [logg for (_, (_, logg, _, _)) in ew_filter_comp]
		FeHs = [FeH for (_, (_, _, FeH, _)) in ew_filter_comp]
		
		scatter!(p_Teff_ew_filter, Ns, abs.(value.(Teffs) .- Teff_ref); yscale=:log10, linewidth=2, ylim=(10^1.6, 10^2.7), xlabel="N (no. iterations)", ylabel="Absolute error / K")
		scatter!(p_logg_ew_filter, Ns, abs.(value.(loggs) .- logg_ref); yscale=:log10, linewidth=2, xlabel="N (no. iterations)", ylabel="Absolute error / cm s^-2")
		scatter!(p_FeH_ew_filter, Ns, abs.(value.(FeHs) .- FeH_ref); yscale=:log10, linewidth=2, xlabel="N (no. iterations)", ylabel="Absolute error / dex")
	end
	
	#p_Teff
	#for (N, (Teff, logg, FeH, α_Fe)) in ew_filter_comp
	#	scatter!([N], [Teff]; label="N="*string(N), markersize=2, markerstrokewidth=0)
	#end
	#p_1_Teff#(k[3], p_1_Teff)
end

# ╔═╡ e3657c40-9eaa-11eb-0e31-a768702858a7
md"Window size = 1"

# ╔═╡ e72eb940-9eaa-11eb-1600-5f4fe79168de
p_Teff_ew_filter

# ╔═╡ ebb9e6b0-9eaa-11eb-3654-bb712597a80f
p_logg_ew_filter

# ╔═╡ ef8baf30-9eaa-11eb-3c4b-b310ec089d2f
p_FeH_ew_filter

# ╔═╡ f36f7910-9eaa-11eb-310a-45ab97855651
md"There seems to be no recognizable pattern. If it were the case that more filter iterations produced better results, I should have seen all the points converging to the smallest error as x increased. That is not observed. Furthermore, the uncertainties are all more or less within range of uncertainty (i.e. step) in grid parameters (2*uncertainty, could swing towards the next step), which means the code can't really produce any more accurate values."

# ╔═╡ f851a020-9eaa-11eb-24ab-156b02dfb2a1
md"N=0 and window size = 1 or window size = 5 should replicate these results. Perhaps window size = 5 helps? Likely not, but let's see."

# ╔═╡ fe207120-9eaa-11eb-08e1-870ace38f2ed
begin
	ew_filter_n0_window_size = Dict()
	
	for filename in filter(x -> endswith(x, "fits"), readdir(sisma_archive_path))[1:20]
		ew_filter_comp = Dict()

		λ, flux, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))

		vrad = opt_data["VRAD"]

		grid_cache = GridCache()

		for window_size=[1,5]
			lines, ew_list = isolate_all_lines_found_in_spectrum(
				vrad,
				keys(full_line_list),
				(λ, flux),
				0,
			    window_size,
				25)

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

			push!(ew_filter_comp, window_size => (Teff, logg, FeH, α_Fe))
		end
		
		push!(ew_filter_n0_window_size, filename => ew_filter_comp)
	end

	ew_filter_n0_window_size
end

# ╔═╡ 01f239a0-9eab-11eb-3c77-29b735566b38
begin
	p_Teff_ew_filter_n0_ws = scatter([], [], label="", legend=nothing)#:outerright)
	p_logg_ew_filter_n0_ws = scatter([], [], label="", legend=nothing)
	p_FeH_ew_filter_n0_ws = scatter([], [], label="", legend=nothing)

	for (filename, ew_filter_comp) in ew_filter_n0_window_size
		_, _, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))

		Teff_ref = opt_data["TEFF"]
		logg_ref = opt_data["LOGG"]
		FeH_ref = opt_data["FEH"]
			
		ew_filter_comp = sort!(collect(ew_filter_comp); by=x -> x[1])
		
		Ns = [N for (N, _) in ew_filter_comp]
		Teffs = [Teff for (_, (Teff, _, _, _)) in ew_filter_comp]
		loggs = [logg for (_, (_, logg, _, _)) in ew_filter_comp]
		FeHs = [FeH for (_, (_, _, FeH, _)) in ew_filter_comp]
		
		scatter!(p_Teff_ew_filter_n0_ws, Ns, abs.(value.(Teffs) .- Teff_ref); yscale=:log10, linewidth=2, ylim=(10^1.6, 10^2.7), xlabel="Window size", ylabel="Absolute error / K")
		scatter!(p_logg_ew_filter_n0_ws, Ns, abs.(value.(loggs) .- logg_ref); yscale=:log10, linewidth=2, xlabel="Window size", ylabel="Absolute error / cm s^-2")
		scatter!(p_FeH_ew_filter_n0_ws, Ns, abs.(value.(FeHs) .- FeH_ref); yscale=:log10, linewidth=2, xlabel="Window size", ylabel="Absolute error / dex")
	end
end

# ╔═╡ 06a47710-9eab-11eb-2ead-4bb75bccb346
p_Teff_ew_filter_n0_ws

# ╔═╡ 0a3e3fa0-9eab-11eb-119b-1f2b6028b29f
p_logg_ew_filter_n0_ws

# ╔═╡ 0d4f0210-9eab-11eb-3fa7-a35ae75a8b7b
p_FeH_ew_filter_n0_ws

# ╔═╡ 1192510e-9eab-11eb-1337-b5932bb7da0e
md"No significant difference."

# ╔═╡ 7bb97b60-9eb3-11eb-13a4-df7c650d2608
begin
	p_Teff_ew_filter_n0_ws_u = scatter([], [], label="", legend=nothing)#:outerright)
	p_logg_ew_filter_n0_ws_u = scatter([], [], label="", legend=nothing)
	p_FeH_ew_filter_n0_ws_u = scatter([], [], label="", legend=nothing)

	for (filename, ew_filter_comp) in ew_filter_n0_window_size
		_, _, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))

		Teff_ref = opt_data["TEFF"]
		logg_ref = opt_data["LOGG"]
		FeH_ref = opt_data["FEH"]
			
		ew_filter_comp = sort!(collect(ew_filter_comp); by=x -> x[1])
		
		Ns = [N for (N, _) in ew_filter_comp]
		Teffs = [Teff for (_, (Teff, _, _, _)) in ew_filter_comp]
		loggs = [logg for (_, (_, logg, _, _)) in ew_filter_comp]
		FeHs = [FeH for (_, (_, _, FeH, _)) in ew_filter_comp]
		
		scatter!(p_Teff_ew_filter_n0_ws_u, Ns, uncertainty.(Teffs); linewidth=2, xlabel="Window size", ylabel="Measurement uncertainty / K")
		scatter!(p_logg_ew_filter_n0_ws_u, Ns, uncertainty.(loggs); linewidth=2, xlabel="Window size", ylabel="Measurement uncertainty / cm s^-2")
		scatter!(p_FeH_ew_filter_n0_ws_u, Ns, uncertainty.(FeHs); linewidth=2, xlabel="Window size", ylabel="Measurement uncertainty / dex")
	end
end

# ╔═╡ de76ed92-9eb4-11eb-2aeb-a9d64be6cceb
p_Teff_ew_filter_n0_ws_u

# ╔═╡ f0d317c0-9eb4-11eb-3d2a-0bf3dcaf4e1d
p_logg_ew_filter_n0_ws_u

# ╔═╡ f8335e30-9eb4-11eb-2abb-b72a2ef0e2c2
p_FeH_ew_filter_n0_ws_u

# ╔═╡ b9aaf4b0-9eb5-11eb-2a08-e1f268953884
md"Uncertainties are not significantly smaller, or they are inexistent, window size = 1 or window size = 5. Perhaps not so with N=0 and N=5."

# ╔═╡ e21861d0-9eb5-11eb-3701-379e2640a589
begin
	p_Teff_ew_filter_n0_5_u = scatter([], [], label="", legend=nothing)#:outerright)
	p_logg_ew_filter_n0_5_u = scatter([], [], label="", legend=nothing)
	p_FeH_ew_filter_n0_5_u = scatter([], [], label="", legend=nothing)

	for (filename, ew_filter_comp) in all_spectra_ew_filter_comp
		_, _, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))

		Teff_ref = opt_data["TEFF"]
		logg_ref = opt_data["LOGG"]
		FeH_ref = opt_data["FEH"]
			
		ew_filter_comp = sort!(collect(ew_filter_comp); by=x -> x[1])
		
		Ns = [N for (N, _) in ew_filter_comp]
		Teffs = [Teff for (_, (Teff, _, _, _)) in ew_filter_comp]
		loggs = [logg for (_, (_, logg, _, _)) in ew_filter_comp]
		FeHs = [FeH for (_, (_, _, FeH, _)) in ew_filter_comp]
		
		plot!(p_Teff_ew_filter_n0_5_u, Ns, uncertainty.(Teffs); linestyle=:dash, linewidth=2, xlabel="No. iterations", ylabel="Teff obtained, uncertainty / K")
		plot!(p_logg_ew_filter_n0_5_u, Ns, uncertainty.(loggs); linestyle=:dash, linewidth=2, xlabel="No. iterations", ylabel="logg obtained, uncertainty / cm s^-2")
		plot!(p_FeH_ew_filter_n0_5_u, Ns, uncertainty.(FeHs); linestyle=:dash, linewidth=2, xlabel="No. iterations", ylabel="[Fe/H] obtained, uncertainty / dex")
	end
end

# ╔═╡ 1406f990-9eb6-11eb-3fe2-b59e1baf1be7
p_Teff_ew_filter_n0_5_u

# ╔═╡ 16f7ff00-9eb6-11eb-3130-b321b24cd955
p_logg_ew_filter_n0_5_u

# ╔═╡ 1adade80-9eb6-11eb-2d94-815be0e86182
p_FeH_ew_filter_n0_5_u

# ╔═╡ Cell order:
# ╠═b2becdce-9eaa-11eb-3e6e-0fb4203c73ec
# ╠═c3c8ae70-9eaa-11eb-2d7d-310fcee9b3ed
# ╠═c7670ae0-9eaa-11eb-3c60-6b9f4b0ce567
# ╠═ca61fb60-9eaa-11eb-2222-dd88d37de4b7
# ╠═ce2d8250-9eaa-11eb-067e-e19e370f8b05
# ╠═d1580e50-9eaa-11eb-1364-55020bd1fe06
# ╠═ddee8f40-9eaa-11eb-23c9-ff7ebe92c868
# ╟─e3657c40-9eaa-11eb-0e31-a768702858a7
# ╠═e72eb940-9eaa-11eb-1600-5f4fe79168de
# ╠═ebb9e6b0-9eaa-11eb-3654-bb712597a80f
# ╠═ef8baf30-9eaa-11eb-3c4b-b310ec089d2f
# ╟─f36f7910-9eaa-11eb-310a-45ab97855651
# ╟─f851a020-9eaa-11eb-24ab-156b02dfb2a1
# ╠═fe207120-9eaa-11eb-08e1-870ace38f2ed
# ╠═01f239a0-9eab-11eb-3c77-29b735566b38
# ╠═06a47710-9eab-11eb-2ead-4bb75bccb346
# ╠═0a3e3fa0-9eab-11eb-119b-1f2b6028b29f
# ╠═0d4f0210-9eab-11eb-3fa7-a35ae75a8b7b
# ╟─1192510e-9eab-11eb-1337-b5932bb7da0e
# ╠═7bb97b60-9eb3-11eb-13a4-df7c650d2608
# ╠═de76ed92-9eb4-11eb-2aeb-a9d64be6cceb
# ╠═f0d317c0-9eb4-11eb-3d2a-0bf3dcaf4e1d
# ╠═f8335e30-9eb4-11eb-2abb-b72a2ef0e2c2
# ╟─b9aaf4b0-9eb5-11eb-2a08-e1f268953884
# ╠═e21861d0-9eb5-11eb-3701-379e2640a589
# ╠═1406f990-9eb6-11eb-3fe2-b59e1baf1be7
# ╠═16f7ff00-9eb6-11eb-3130-b321b24cd955
# ╠═1adade80-9eb6-11eb-2d94-815be0e86182
