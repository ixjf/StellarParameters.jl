### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 98072c10-9fa1-11eb-0da6-5bc25317b68a
begin
	using Pkg
	Pkg.activate("..")
	
	using Plots
	gr();
	
	using Measurements: value, uncertainty, Measurement, ±
	
	using FFTW, DSP
	
	using Statistics
	
	using Format
	
	using Interpolations
end

# ╔═╡ de00af50-ee65-4280-b2e9-8d448f5246c6
include("../src/StellarParameters.jl");

# ╔═╡ 959dceb1-7237-448a-ae1d-76cb85dca9c3
full_line_list, ew_list_Sun = read_line_list()

# ╔═╡ 6394896f-83fc-47f4-95ba-72c158888829
multiplet_list = group_lines_into_multiplets(full_line_list)

# ╔═╡ 7399902e-0f5d-4347-b4a1-188a96fcfed8
spectra_archive_path = joinpath(@__DIR__, "../../data/estrelas_a_analisar/")

# ╔═╡ 3aa0fdcf-2627-46a9-9d49-c86514f4ec27
md"# Determining stellar parameters"

# ╔═╡ ca6d2c9e-e0cb-4a40-9d7c-104ca912a98c
begin
	results = []
	
	grid_cache = GridCache()
	
	for filename in ["estrela1.fits", "estrela2_vcor.fits"]
		λ, flux, Δx = read_spectrum_unknown_source_prof(joinpath(spectra_archive_path, filename))
		 
		lines, ew_list = isolate_all_lines_found_in_spectrum(
			0.0,
			keys(full_line_list),
			(λ, flux),
			5, # already determined in 2_lines.jl
			5, # already determined in 2_lines.jl
			12.5) # already determined in 2_lines.jl

		best_combo_i, best_combo_j, _ = find_best_multiplet_combo_for_Texc_estimate(
			weighted_median, 
			multiplet_list, 
			Dict([λ_c_approx => (λ_c_approx,) for (λ_c_approx,_) in lines]), 
			ew_list_Sun,
			2)

		# Spectra are lower res, harder to find lines. can't use combo (1,4)
		_, _, Texc = find_best_multiplet_combo_for_Texc_estimate(
			weighted_median, 
			[multiplet_list[best_combo_i], multiplet_list[best_combo_j]],
			lines,
			ew_list,
			2)

		params = match_ew_against_grid(ew_list, Texc, grid_cache)
		
		isnothing(params) && continue

		Teff, logg, FeH, α_Fe = params

		push!(
			results, 
			filename => (λ, 
				flux, 
				lines, 
				ew_list, 
		      	best_combo_i,
		      	best_combo_j,
				Texc, 
				Teff, 
				logg, 
				FeH, 
				α_Fe))
	end
end

# ╔═╡ 60b604c6-15e1-44f6-a369-2f31dea63dfa
md"While it seemed best to use combo (1,4) before, we find that for these lower resolution spectra, some multiplets don't have enough lines well defined, so we do have to use the best combo for the Sun."

# ╔═╡ db2ae75b-cd61-40ec-ba80-811851ad5b5f
md"**estrela1.fits**

Best multiplet combo: $(results[1][2][5], results[1][2][6])

Texc: $(value(results[1][2][7])) ± $(uncertainty(results[1][2][7]))

Teff: $(value(results[1][2][8])) ± $(uncertainty(results[1][2][8]))

logg: $(value(results[1][2][9])) ± $(uncertainty(results[1][2][9]))

[Fe/H]: $(value(results[1][2][10])) ± $(uncertainty(results[1][2][10]))

[α/Fe]: $(value(results[1][2][11])) ± $(uncertainty(results[1][2][11]))"

# ╔═╡ dce9b2d7-6164-43c3-a2ad-c8472e74bb6a
md"**estrela2_vcor.fits**

Best multiplet combo: $(results[2][2][5], results[2][2][6])

Texc: $(value(results[2][2][7])) ± $(uncertainty(results[2][2][7]))

Teff: $(value(results[2][2][8])) ± $(uncertainty(results[2][2][8]))

logg: $(value(results[2][2][9])) ± $(uncertainty(results[2][2][9]))

[Fe/H]: $(value(results[2][2][10])) ± $(uncertainty(results[2][2][10]))

[α/Fe]: $(value(results[2][2][11])) ± $(uncertainty(results[2][2][11]))"

# ╔═╡ 62872c0c-1b41-4908-9d80-5d8cc8e945b3
function fft_line_log10(λ, flux, continuum, N)
	dλ = λ[2] - λ[1]
	
	wave = [repeat([λ[begin]], N÷2), λ..., repeat([λ[end]], N÷2)]
	flux = [repeat([flux[begin]], N÷2)..., flux..., repeat([flux[end]], N÷2)...]
	
	ft = fft(flux .- continuum)
	ps = log10.(abs.(ft))
	freq = DSP.fftfreq(length(flux), 1/dλ)
	
	ps = ps[freq .> 0]
	freq = freq[freq .> 0]
	
	(freq, ps)
end

# ╔═╡ 4b8403f7-28b6-430a-9763-d7172f6f4839
function line_data(spectrum_i, line_i)
	k = collect(keys(results[spectrum_i][2][3]))
	line = (results[spectrum_i][2][3][k[line_i]])	
end

# ╔═╡ e522a12b-c830-499b-9a78-c0ef69351506
function plot_line(line)
	λ = line[2]
	flux = line[3]
	scatter(λ, flux; legend=nothing)
end

# ╔═╡ f5af9039-71fd-43da-b5af-4a390031d65b
function plot_fft_line(line, continuum, N, noise_level, xlim=(-Inf, +Inf); σ₁lim=missing)
	p = scatter([], []; label="")
	
	λ_c = line[1]
	λ = line[2]
	flux = line[3]
	
	freq, ps = fft_line_log10(λ, flux, continuum, N)
	
	plot!(p, freq, ps; markersize=1, label="FFT (λ_c: $(value(λ_c)))", xlabel="Frequência / Hz", ylabel="log10(fluxo) / dex", xlim=xlim)
	hline!(p, [noise_level]; label="Nível do ruído")

	if !ismissing(σ₁lim)
		I = σ₁lim[1] .<= freq .<= σ₁lim[2]
		σ₁_y, σ₁_arridx = findmin(ps[I])
		σ₁ = freq[I][σ₁_arridx]
		scatter!(p, [σ₁], [σ₁_y]; label="σ₁ = $(σ₁ ± (freq[2] - freq[1]))", markersize=2)
	end
	
	p
end

# ╔═╡ 1d74d647-2e9c-4381-b265-47d7c309ae1d
Δλₘ(σ₁) = 0.660/σ₁

# ╔═╡ ada27e13-a07e-4212-8cde-e6dd1494c916
vsini(σ₁, λ_c) = (Δλₘ(σ₁)*c/λ_c)*1e-3 # km/s

# ╔═╡ 72327cd9-01bc-4f14-a8e6-8e348d6aae21
md"# vsini: **estrela1.fits**"

# ╔═╡ 45203e4a-aa5f-4d42-ab79-d8dc8231c434
md"Continuum, noise level, and σₙ are obtained visually. Only used lines with Gaussian shape and approx constant continuum level (in time domain) & constant noise level (in frequency domain). Large peaks below the estimated noise level are assumed to be flukes. The zero should be very, very close to the noise level, or I discard the line."

# ╔═╡ 1f49530c-7650-4ec8-8ac4-5105cf793b8f
σ₁_estrela1 = Tuple{Float64, Measurement{Float64}}[]

# ╔═╡ ff627fec-82af-4e34-a86a-ffc6f9128c7d
plot_fft_line(line_data(1, 1), 1.5, 1000, -2.4; σ₁lim=(13.9, 14.2))

# ╔═╡ d3a3446e-c375-416c-aefa-61b130e4de29
push!(σ₁_estrela1, (5618.645639486446, 14.047 ± 0.093))

# ╔═╡ d867e55f-b332-4fb2-a1eb-40e861c4cecd
plot_fft_line(line_data(1, 10), 1.05, 1000, -3.05; σ₁lim=(30.05, 30.06))

# ╔═╡ 044ae049-735d-4eca-9b31-d159156e6387
push!(σ₁_estrela1, (6858.165240478099, 30.058 ± 0.095))

# ╔═╡ 41a35fa7-3665-4f79-9b92-42ac9da54805
plot_fft_line(line_data(1, 18), 1.36, 1000, -2.71; σ₁lim=(15.94, 15.97))

# ╔═╡ 8f6cf400-f5bf-40e8-b188-cd303eef367d
push!(σ₁_estrela1, (6710.340049661337, 15.952 ± 0.093))

# ╔═╡ 985b995b-0779-4b8b-96c7-8e8b8a540592
plot_fft_line(line_data(1, 26), 1.40, 1000, -2.48; σ₁lim=(10.65, 10.66))

# ╔═╡ ec25e71c-2637-496b-b1b9-c60c5a71dc40
push!(σ₁_estrela1, (5927.810086326589, 10.655 ± 0.093))

# ╔═╡ 467ed682-b96e-49a2-8ad7-20fdeac7a7ec
plot_fft_line(line_data(1, 33), 1.36, 1000, -2.8; σ₁lim=(14.19, 14.195))

# ╔═╡ e387bb75-9fbe-49d0-b3e1-edee94a9cb97
push!(σ₁_estrela1, (5815.232975922436, 14.191 ± 0.096))

# ╔═╡ a2bf80b9-b60d-4edd-87ec-b68a5acfb8df
begin
	# https://discourse.julialang.org/t/how-to-draw-a-rectangular-region-with-plots-jl/1719/8
	rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
	
	p_vsini_estrela_1 = plot_fft_line(line_data(1, 54), 1.32, 1000, -2.6; σ₁lim=(17.298, 17.302))
	plot!(p_vsini_estrela_1, rectangle(10,0.7,40,-2.95), opacity=.2; fillcolor=:grey, label="Zona do nível do ruído")
	savefig(p_vsini_estrela_1, "plot_vsini_estrela_1.pdf")
	p_vsini_estrela_1
end

# ╔═╡ 9430188d-1178-4411-9e9f-07e7100b28e3
push!(σ₁_estrela1, (5987.080608350745, 17.3 ± 0.093))

# ╔═╡ 4ee26d24-f917-4948-913d-88a0ba8f7c51
plot_fft_line(line_data(1, 57), 1.315, 1000, -2.32; σ₁lim=(22.18, 22.19))

# ╔═╡ 1e453913-5fbf-4d2c-a239-b6eb75ace2d5
push!(σ₁_estrela1, (6024.0744990292305, 22.182 ± 0.092))

# ╔═╡ 3474b583-58fc-4f96-b2bd-03e0ecb61008
plot_fft_line(line_data(1, 67), 1.30, 1000, -2.51; σ₁lim=(6.717, 6.718))

# ╔═╡ a35eb951-ef85-40d0-936d-ebb1f1321ee5
push!(σ₁_estrela1, (6096.679536240571, 6.718 ± 0.095))

# ╔═╡ 85a58f36-09c5-4d7b-bfcd-6f3e9ab8a6a6
vsini_mean_estrela1 = mean(map(x -> vsini(x[2], x[1]), σ₁_estrela1))

# ╔═╡ f0a97221-9c75-437b-baea-9f7b80567a29
md"# vsini: estrela2_vcor.fits"

# ╔═╡ 24ea1dc5-25c4-4760-9fa4-a31d8c5fd207
σ₁_estrela2 = Tuple{Float64, Measurement{Float64}}[]

# ╔═╡ 972bdda9-4371-4fe1-89b7-86f0bb5ed163
begin
	#plot_line(line_data(2, 1))
	#hline!([5150])
	#push!(σ1_estrela1, (6024.0744990292305, 22.186 ± 0.002))
	plot_fft_line(line_data(2, 7), 5150, 1000, 1.75)#.32, (22.185, 22.187)) # σ₁ ≈ 22.186
end

# fuck - no noise level?! but we do know the first zero is around the first minimum, right?

# ╔═╡ 2017e0a2-cd9b-4a1a-8106-52b032b022c4
md"Unable to determine vsini for estrela2_vcor.fits"

# ╔═╡ d98b8d9a-6ab1-4d9a-9ccc-77e3e63103d4
md"# Comparison of observed vs best-fit synthetic spectra"

# ╔═╡ 490f8232-6959-4c4f-a44a-531da801f63f
flux_ambre_folder = "GES_UVESRed580_deltaAlphaFe+0.0_fits"

# ╔═╡ b0b08940-50da-48bd-a656-d0c82101d358
ambre_proj_path = "../../data/ambre_proj/"

# ╔═╡ c53b6180-4348-4148-b72e-52841d2b4abe
# FIXME: hack, I shouldn't manually insert path to grid spectrum here.
wave = read_spectrum_ambre_grid(joinpath(ambre_proj_path, "GES_UVESRed580_Lambda.fits"))

# ╔═╡ a2c356b2-a249-4da9-a229-fcb09df28c4d
flux_star_file(Teff, logg, z, a) = "p$(Teff)_g$(logg)_m0.0_t01_z$(z)_a$(a).GES4750.fits"

# ╔═╡ 1d15c6cc-fd62-4824-9c87-fa5535e709e6
begin
	R_obs = 60000
	ϵ = EPSILON
end

# ╔═╡ c999b35f-55ac-459d-af19-27c932eca48d
md"## estrela1.fits"

# ╔═╡ c44f5195-0bfa-4d44-9c21-f1dad18d7e6f
lower_bound_1, upper_bound_1 = 5700, 5730 # Å

# ╔═╡ 802e490f-76b5-4852-a2cc-951fe69a87b3
λ_estrela1_full, flux_estrela1_full, Δx_estrela1 = read_spectrum_unknown_source_prof(joinpath(spectra_archive_path, "estrela1.fits"))

# ╔═╡ 68ba0c3a-f9de-4260-a141-5ef93eb89ccb
begin
	Teff_estrela1_str = Int64(value(results[1][2][8]))
	logg_estrela1_str = cfmt("%+.1f", value(results[1][2][9]))
	z_estrela1_str = cfmt("%+.2f", value(results[1][2][10]))
	a_estrela1_str = cfmt("%+.2f", value(results[1][2][11]))

	flux_estrela1_bestfit_full = read_spectrum_ambre_grid(joinpath(ambre_proj_path, flux_ambre_folder, flux_star_file(Teff_estrela1_str, logg_estrela1_str, z_estrela1_str, a_estrela1_str)))
end

# ╔═╡ d8836d7a-37d7-46a9-9c97-15090e04780a
md"### Taking interval that we want to visualize in spectrum"

# ╔═╡ 0ccb2ad3-b69b-440f-95b2-5730902ab6a0
begin
	wave_estrela1_bestfit = wave[@. lower_bound_1 <= wave <= upper_bound_1]
	flux_estrela1_bestfit = flux_estrela1_bestfit_full[@. lower_bound_1 <= wave <= upper_bound_1]
	λ_estrela1 = λ_estrela1_full[@. lower_bound_1 <= λ_estrela1_full <= upper_bound_1]
	flux_estrela1 = flux_estrela1_full[@. lower_bound_1 <= λ_estrela1_full <= upper_bound_1]
end

# ╔═╡ be15f58e-d743-49a9-8596-474e40875b64
md"### Interpolating grid spectrum to same Δλ as obs spectrum"

# ╔═╡ c42dbd1c-239c-4c10-a90e-426ba0c1323b
begin
	itp = LinearInterpolation(wave_estrela1_bestfit, flux_estrela1_bestfit)
	flux_estrela1_bestfit_itp = itp.(λ_estrela1)
	wave_estrela1_bestfit_itp = λ_estrela1
end

# ╔═╡ 2187b1ac-11be-4383-900c-0a947d15f95c
md"### Normalizing observed spectrum to continuum level"

# ╔═╡ e4776dd5-361b-49b7-8013-25e6445587d8
begin
	p_estrela1_normalizacao = plot(λ_estrela1, flux_estrela1; markersize=1, label="Espetro real", legend=:bottomright, xlim=(lower_bound_1, upper_bound_1))
	y0, y1 = 1.38, 1.385
	plot!(p_estrela1_normalizacao, [lower_bound_1, upper_bound_1], [y0, y1]; label="Nível do contínuo")
	xlabel!(p_estrela1_normalizacao, "Comprimento de onda / Å")
	ylabel!(p_estrela1_normalizacao, "Fluxo / W/m²")
	savefig(p_estrela1_normalizacao, "estrela1_normalizacao.pdf")
	p_estrela1_normalizacao
end

# ╔═╡ c2546bed-788e-4241-a4fa-bb5c82fff028
begin
	m = (y1 - y0)/(upper_bound_1 - lower_bound_1)
	b = y0 - m*lower_bound_1
	
	continuum_estrela1 = [m*x + b for x in λ_estrela1]
	
	flux_estrela1_normalized = flux_estrela1 ./ continuum_estrela1
	plot(λ_estrela1, flux_estrela1_normalized; legend=nothing)
end

# ╔═╡ 00a1a885-fe9c-4e9c-a609-772220416a14
md"### Convolving with rotational + instrumental profiles"

# ╔═╡ 5ef8990d-ba2f-4890-84cd-fbf10f17abbc
begin
	wave_estrela1_bestfit_sim, flux_estrela1_bestfit_sim = convolve_simulated_profile(R_obs, vsini_mean_estrela1, ϵ, mean([upper_bound_1, lower_bound_1]), repeat([1.0], length(wave_estrela1_bestfit_itp)), wave_estrela1_bestfit_itp, flux_estrela1_bestfit_itp)
end

# ╔═╡ 0fc4f855-711b-4c97-9e29-92cb907ebd5b
begin
	plot(wave_estrela1_bestfit_sim, flux_estrela1_bestfit_sim; label="Simulated")
	plot!(wave_estrela1_bestfit_itp, flux_estrela1_bestfit_itp; label="Raw", legend=:bottom)
end

# ╔═╡ 1513ef26-7844-4cca-9369-3f0538d8d540
md"### Final result"

# ╔═╡ 26ad8858-ec8b-4bdf-bad4-babf78ee3c46
begin
	p_final_result_estrela_1 = plot(wave_estrela1_bestfit_sim, flux_estrela1_bestfit_sim; label="Melhor ajuste", legend=:bottomright)
	plot!(p_final_result_estrela_1, λ_estrela1, flux_estrela1_normalized; label="Espetro real")
	xlabel!(p_final_result_estrela_1, "Comprimento de onda / Å")
	ylabel!(p_final_result_estrela_1, "Fluxo / W/m²")
	savefig(p_final_result_estrela_1, "resultado_final_estrela_1.pdf")
	p_final_result_estrela_1
end

# ╔═╡ ac1012b7-bc23-4ad4-ba39-e6c87db9494a
md"The difference in signal intensity is at least partially due to normalization, it kinda fucks that up.

TODO: Why do some lines not fit?!"

# ╔═╡ 6c1fd78d-5973-4f35-9999-a1608c885ddc
md"## estrela2_vcor.fits"

# ╔═╡ 2f0d7679-2cfd-4378-b4b6-d0396ae75c22
lower_bound_2, upper_bound_2 = 6300, 6330 # Å

# ╔═╡ ddea0b06-d28f-4c36-9dcb-f30a96449f0f
λ_estrela2_full, flux_estrela2_full, Δx_estrela2 = read_spectrum_unknown_source_prof(joinpath(spectra_archive_path, "estrela2_vcor.fits"))

# ╔═╡ eb7a7c9a-bc07-4bc3-afef-e5c3029b38e5
begin
	Teff_estrela2_str = Int64(value(results[2][2][8]))
	logg_estrela2_str = cfmt("%+.1f", value(results[2][2][9]))
	z_estrela2_str = cfmt("%+.2f", value(results[2][2][10]))
	a_estrela2_str = cfmt("%+.2f", value(results[2][2][11]))

	flux_estrela2_bestfit_full = read_spectrum_ambre_grid(joinpath(ambre_proj_path, flux_ambre_folder, flux_star_file(Teff_estrela2_str, logg_estrela2_str, z_estrela2_str, a_estrela2_str)))
end

# ╔═╡ 71a00e02-341e-41e0-845f-8213c2d694cb
md"### Taking interval that we want to visualize in spectrum"

# ╔═╡ 35f97c8c-d77a-4fe3-b12a-20e887b7c485
begin
	wave_estrela2_bestfit = wave[@. lower_bound_2 <= wave <= upper_bound_2]
	flux_estrela2_bestfit = flux_estrela2_bestfit_full[@. lower_bound_2 <= wave <= upper_bound_2]
	λ_estrela2 = λ_estrela2_full[@. lower_bound_2 <= λ_estrela2_full <= upper_bound_2]
	flux_estrela2 = flux_estrela2_full[@. lower_bound_2 <= λ_estrela2_full <= upper_bound_2]
end

# ╔═╡ f4dbd17f-9587-4ed6-8d24-afa2a14a2e9c
md"### Interpolating grid spectrum to same Δλ as obs spectrum"

# ╔═╡ d9a3af6e-48fa-46f1-81b0-b268db043657
begin
	itp2 = LinearInterpolation(wave_estrela2_bestfit, flux_estrela2_bestfit)
	flux_estrela2_bestfit_itp = itp2.(λ_estrela2)
	wave_estrela2_bestfit_itp = λ_estrela2
end

# ╔═╡ 09368cfe-89a4-4dac-bff9-58ef30ac469d
md"### Normalizing observed spectrum to continuum level"

# ╔═╡ 90201c5b-78c8-47c2-ad2c-f29926c4d1ba
begin
	p_estrela2_normalizacao = plot(λ_estrela2, flux_estrela2; markersize=1, label="Espetro real", legend=:bottom, xlim=(lower_bound_2, upper_bound_2))
	y0_2, y1_2 = 5000, 5350
	plot!(p_estrela2_normalizacao, [lower_bound_2, upper_bound_2], [y0_2, y1_2]; label="Nível do contínuo")
	xlabel!(p_estrela2_normalizacao, "Comprimento de onda / Å")
	ylabel!(p_estrela2_normalizacao, "Fluxo / W/m²")
	savefig(p_estrela2_normalizacao, "estrela2_normalizacao.pdf")
	p_estrela2_normalizacao
end

# ╔═╡ 21cd87b3-3dd1-4bbc-a1f8-b673b9d75bf5
begin
	m_2 = (y1_2 - y0_2)/(upper_bound_2 - lower_bound_2)
	b_2 = y0_2 - m_2*lower_bound_2
	
	continuum_estrela2 = [m_2*x + b_2 for x in λ_estrela2]
	
	flux_estrela2_normalized = flux_estrela2 ./ continuum_estrela2
	plot(λ_estrela2, flux_estrela2_normalized; legend=nothing)
end

# ╔═╡ 47b85320-d9a9-4df4-877b-d26f1798d21e
md"### Convolving with rotational + instrumental profiles"

# ╔═╡ 965dee48-f698-47bc-83fe-4d84c9aa7293
begin
	wave_estrela2_bestfit_sim, flux_estrela2_bestfit_sim = convolve_simulated_profile(R_obs, missing, ϵ, mean([upper_bound_2, lower_bound_2]), repeat([1.0], length(wave_estrela2_bestfit_itp)), wave_estrela2_bestfit_itp, flux_estrela2_bestfit_itp)
end

# ╔═╡ 867588d1-b18c-42a4-957e-8797ec678af0
begin
	plot(wave_estrela2_bestfit_sim, flux_estrela2_bestfit_sim; label="Simulated")
	plot!(wave_estrela2_bestfit_itp, flux_estrela2_bestfit_itp; label="Raw", legend=:bottom)
end

# ╔═╡ d3d5474e-2b73-49e6-a07f-1f49fe1c2205
begin
	p_final_result_estrela_2 = plot(wave_estrela2_bestfit_sim, flux_estrela2_bestfit_sim; label="Melhor ajuste", legend=:bottomright)
	plot!(p_final_result_estrela_2, λ_estrela2, flux_estrela2_normalized; label="Espetro real")
	xlabel!(p_final_result_estrela_2, "Comprimento de onda / Å")
	ylabel!(p_final_result_estrela_2, "Fluxo / W/m²")
	savefig(p_final_result_estrela_2, "resultado_final_estrela_2.pdf")
	p_final_result_estrela_2
end

# ╔═╡ Cell order:
# ╠═98072c10-9fa1-11eb-0da6-5bc25317b68a
# ╠═de00af50-ee65-4280-b2e9-8d448f5246c6
# ╠═959dceb1-7237-448a-ae1d-76cb85dca9c3
# ╠═6394896f-83fc-47f4-95ba-72c158888829
# ╠═7399902e-0f5d-4347-b4a1-188a96fcfed8
# ╟─3aa0fdcf-2627-46a9-9d49-c86514f4ec27
# ╠═ca6d2c9e-e0cb-4a40-9d7c-104ca912a98c
# ╟─60b604c6-15e1-44f6-a369-2f31dea63dfa
# ╠═db2ae75b-cd61-40ec-ba80-811851ad5b5f
# ╠═dce9b2d7-6164-43c3-a2ad-c8472e74bb6a
# ╠═62872c0c-1b41-4908-9d80-5d8cc8e945b3
# ╠═4b8403f7-28b6-430a-9763-d7172f6f4839
# ╠═e522a12b-c830-499b-9a78-c0ef69351506
# ╠═f5af9039-71fd-43da-b5af-4a390031d65b
# ╠═1d74d647-2e9c-4381-b265-47d7c309ae1d
# ╠═ada27e13-a07e-4212-8cde-e6dd1494c916
# ╟─72327cd9-01bc-4f14-a8e6-8e348d6aae21
# ╟─45203e4a-aa5f-4d42-ab79-d8dc8231c434
# ╠═1f49530c-7650-4ec8-8ac4-5105cf793b8f
# ╠═ff627fec-82af-4e34-a86a-ffc6f9128c7d
# ╠═d3a3446e-c375-416c-aefa-61b130e4de29
# ╠═d867e55f-b332-4fb2-a1eb-40e861c4cecd
# ╠═044ae049-735d-4eca-9b31-d159156e6387
# ╠═41a35fa7-3665-4f79-9b92-42ac9da54805
# ╠═8f6cf400-f5bf-40e8-b188-cd303eef367d
# ╠═985b995b-0779-4b8b-96c7-8e8b8a540592
# ╠═ec25e71c-2637-496b-b1b9-c60c5a71dc40
# ╠═467ed682-b96e-49a2-8ad7-20fdeac7a7ec
# ╠═e387bb75-9fbe-49d0-b3e1-edee94a9cb97
# ╠═a2bf80b9-b60d-4edd-87ec-b68a5acfb8df
# ╠═9430188d-1178-4411-9e9f-07e7100b28e3
# ╠═4ee26d24-f917-4948-913d-88a0ba8f7c51
# ╠═1e453913-5fbf-4d2c-a239-b6eb75ace2d5
# ╠═3474b583-58fc-4f96-b2bd-03e0ecb61008
# ╠═a35eb951-ef85-40d0-936d-ebb1f1321ee5
# ╠═85a58f36-09c5-4d7b-bfcd-6f3e9ab8a6a6
# ╟─f0a97221-9c75-437b-baea-9f7b80567a29
# ╠═24ea1dc5-25c4-4760-9fa4-a31d8c5fd207
# ╠═972bdda9-4371-4fe1-89b7-86f0bb5ed163
# ╟─2017e0a2-cd9b-4a1a-8106-52b032b022c4
# ╟─d98b8d9a-6ab1-4d9a-9ccc-77e3e63103d4
# ╠═490f8232-6959-4c4f-a44a-531da801f63f
# ╠═b0b08940-50da-48bd-a656-d0c82101d358
# ╠═c53b6180-4348-4148-b72e-52841d2b4abe
# ╠═a2c356b2-a249-4da9-a229-fcb09df28c4d
# ╠═1d15c6cc-fd62-4824-9c87-fa5535e709e6
# ╟─c999b35f-55ac-459d-af19-27c932eca48d
# ╠═c44f5195-0bfa-4d44-9c21-f1dad18d7e6f
# ╠═802e490f-76b5-4852-a2cc-951fe69a87b3
# ╠═68ba0c3a-f9de-4260-a141-5ef93eb89ccb
# ╟─d8836d7a-37d7-46a9-9c97-15090e04780a
# ╠═0ccb2ad3-b69b-440f-95b2-5730902ab6a0
# ╟─be15f58e-d743-49a9-8596-474e40875b64
# ╠═c42dbd1c-239c-4c10-a90e-426ba0c1323b
# ╟─2187b1ac-11be-4383-900c-0a947d15f95c
# ╠═e4776dd5-361b-49b7-8013-25e6445587d8
# ╠═c2546bed-788e-4241-a4fa-bb5c82fff028
# ╟─00a1a885-fe9c-4e9c-a609-772220416a14
# ╠═5ef8990d-ba2f-4890-84cd-fbf10f17abbc
# ╠═0fc4f855-711b-4c97-9e29-92cb907ebd5b
# ╟─1513ef26-7844-4cca-9369-3f0538d8d540
# ╠═26ad8858-ec8b-4bdf-bad4-babf78ee3c46
# ╟─ac1012b7-bc23-4ad4-ba39-e6c87db9494a
# ╟─6c1fd78d-5973-4f35-9999-a1608c885ddc
# ╠═2f0d7679-2cfd-4378-b4b6-d0396ae75c22
# ╠═ddea0b06-d28f-4c36-9dcb-f30a96449f0f
# ╠═eb7a7c9a-bc07-4bc3-afef-e5c3029b38e5
# ╟─71a00e02-341e-41e0-845f-8213c2d694cb
# ╠═35f97c8c-d77a-4fe3-b12a-20e887b7c485
# ╟─f4dbd17f-9587-4ed6-8d24-afa2a14a2e9c
# ╠═d9a3af6e-48fa-46f1-81b0-b268db043657
# ╟─09368cfe-89a4-4dac-bff9-58ef30ac469d
# ╠═90201c5b-78c8-47c2-ad2c-f29926c4d1ba
# ╠═21cd87b3-3dd1-4bbc-a1f8-b673b9d75bf5
# ╟─47b85320-d9a9-4df4-877b-d26f1798d21e
# ╠═965dee48-f698-47bc-83fe-4d84c9aa7293
# ╠═867588d1-b18c-42a4-957e-8797ec678af0
# ╠═d3d5474e-2b73-49e6-a07f-1f49fe1c2205
