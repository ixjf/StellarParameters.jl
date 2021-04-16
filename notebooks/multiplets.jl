### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 40bbc780-8f17-11eb-343e-e1c079a00e28
begin
	using Pkg
	Pkg.activate("..")
	
	using Plots, StatsPlots
	gr()
	
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
			keys(full_line_list),
			(λ, flux))
		
		push!(d, filename => (λ, flux, opt_data, lines, ew_list))
	end
		
	d
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
end

# ╔═╡ 5cb65bae-8ff0-11eb-1d7f-79ea6ef8cf6b
begin
	p = scatter([], []; label="")
	
	scatter!(Teffs, map(x -> x[3], unweighted_multiplet_data); label="Unweighted median", markercolor=:grey22, markersize=2.5, markerstrokecolor=:grey22)
	scatter!(Teffs, map(x -> x[3], weighted_multiplet_data); label="Weighted median", markercolor=:royalblue2, markersize=2.5, markerstrokecolor=:royalblue2)
	
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

# ╔═╡ 353cfb90-9d22-11eb-25b1-fdb1e6bf014c
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
				N)

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

# ╔═╡ a10a4380-9d24-11eb-3215-11828ecd1642
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

# ╔═╡ 60576220-9e94-11eb-24c0-ed61e8eeaf25
md"Window size = 1"

# ╔═╡ 3c843260-9d3b-11eb-060b-e3a7fcc36502
p_Teff_ew_filter

# ╔═╡ 3e634b40-9d3e-11eb-106e-2760d026e94e
p_logg_ew_filter

# ╔═╡ 460814c0-9d3e-11eb-0eee-d12a0da82fc3
p_FeH_ew_filter

# ╔═╡ 92a1e9f0-9d3e-11eb-095e-b32071933ca5
md"There seems to be no recognizable pattern. If it were the case that more filter iterations produced better results, I should have seen all the points converging to the smallest error as x increased. That is not observed. Furthermore, the uncertainties are all more or less within range of uncertainty (i.e. step) in grid parameters (2*uncertainty, could swing towards the next step), which means the code can't really produce any more accurate values."

# ╔═╡ e76f6dc0-9e94-11eb-02f7-83caaf29b4d6
md"N=0 and window size = 1 or window size = 5 should replicate these results. Perhaps window size = 5 helps? Likely not, but let's see."

# ╔═╡ 16c706ee-9e95-11eb-18aa-57a9379c3ff4
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
			    window_size)

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

# ╔═╡ 4862c3c0-9e95-11eb-3e16-7d37c85c5935
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

# ╔═╡ 6b4b8840-9e95-11eb-0d43-cfeddf3508b4
p_Teff_ew_filter_n0_ws

# ╔═╡ 6eaeb070-9e95-11eb-133f-f59df21a7e62
p_logg_ew_filter_n0_ws

# ╔═╡ 7418f340-9e95-11eb-030f-39e81ddcccc2
p_FeH_ew_filter_n0_ws

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
# ╠═353cfb90-9d22-11eb-25b1-fdb1e6bf014c
# ╠═a10a4380-9d24-11eb-3215-11828ecd1642
# ╟─60576220-9e94-11eb-24c0-ed61e8eeaf25
# ╠═3c843260-9d3b-11eb-060b-e3a7fcc36502
# ╠═3e634b40-9d3e-11eb-106e-2760d026e94e
# ╠═460814c0-9d3e-11eb-0eee-d12a0da82fc3
# ╟─92a1e9f0-9d3e-11eb-095e-b32071933ca5
# ╟─e76f6dc0-9e94-11eb-02f7-83caaf29b4d6
# ╠═16c706ee-9e95-11eb-18aa-57a9379c3ff4
# ╠═4862c3c0-9e95-11eb-3e16-7d37c85c5935
# ╠═6b4b8840-9e95-11eb-0d43-cfeddf3508b4
# ╠═6eaeb070-9e95-11eb-133f-f59df21a7e62
# ╠═7418f340-9e95-11eb-030f-39e81ddcccc2
