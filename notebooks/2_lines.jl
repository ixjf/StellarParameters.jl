### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 0ba701f0-9006-11eb-08c9-6183e18fa170
begin
	using Pkg
	Pkg.activate("..")
	
	using Plots
	gr()
	
	using Measurements: value, uncertainty, Measurement
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
	# https://stackoverflow.com/a/58140365
# NOTE: This isn't really necessary except for plotting luminosity as a function of ξ
function cumtrapz(X::T, Y::T) where {T <: AbstractVector}
    # Check matching vector length
    # @assert length(X) == length(Y)
    # Initialize Output
    out = similar(X)
    out[1] = 0
    # Iterate over arrays
    for i in 2:length(X)
        out[i] = out[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
    end
    # Return output
    out
end
	
	ew_list = []

	for filename in filter(x -> endswith(x, "fits"), readdir(sisma_archive_path))
		λ, flux, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))
	
		vrad = opt_data["VRAD"]
	
		lines, ew_list_spec = isolate_all_lines_found_in_spectrum(
			vrad,
			keys(full_line_list),
			(λ, flux),
	    	5,
	    	5,
			+Inf)
		
		for (k,v) in ew_list_spec
			λ_line, flux_line = lines[k][2], lines[k][3]
			
			fit_result = fit_line_as_gaussian((λ_line, flux_line), value(lines[k][1]), +Inf)

			isnothing(fit_result) && continue # ???
			
        	λ_c_exact, σ, A, B, W_λ = fit_result
			
			ew_area = -(cumtrapz(λ_line, flux_line .- value(B))[end])
			push!(ew_list, k => (v, ew_area))
		end
	end
	
	ew_list
end

# ╔═╡ 5d4d630e-0a68-415b-a9cd-8650c6875914
scatter(map(x -> x[1], ew_list), [uncertainty(W_λ)/value(W_λ)*100 for (λ_c_approx, (W_λ, _)) in ew_list]; label="", xlabel="λ_c / Å", ylabel="W_λ % relative uncertainty", ylim=(0, 100))

# ╔═╡ c711a7d4-af33-4445-823b-514c1dbcbf7b
begin
	ew_residuals = [W_λ - W_λ_area for (k, (W_λ, W_λ_area)) in ew_list]
	scatter(map(x -> x[1], ew_list), value.(ew_residuals); label="", xlabel="λ_c / Å", ylabel="W_λ - W_λ_area residuals", ylim=(-0.025, 0.025))
end

# ╔═╡ f5a8aced-e91e-4b63-9eeb-8f7d32d5dc9b
md"The graph above combines the lines found in all ~100 SISMA spectra WITH % unc <= 100%.

The majority of the lines have an EW with % relative uncertainty < 12.5%. This suggests that, in general, up to ~ 12.5% uncertainty is expected. We can assume that any uncertainties higher than that are flukes and thus disregard them."

# ╔═╡ 076fdfdf-096d-4b66-ae81-081d5864c1ec
md"By playing around with the number of filter iterations and window size, we also notice that:

- N = 0, window size = 1 (no filtering at all) => really high uncertainties in many cases
- N = 0, window size = 5 => most values are < 25%
- N > 0, window size = 5 => diminishing returns as N increases
- N > 0, window size = 1 => almost as bad as N = 0, window size = 1

So we can say that window size matters more than number of iterations, and is also more efficient. But N = 5, window size = 5 is only 0.6s slower (total time for all spectra), so it's worth it."

# ╔═╡ 5a3553a3-4af0-4ef3-b8b6-a7e635184636
md"TODO: plot graph of increasing number of lines detected with unc < 12.5% as N increases"

# ╔═╡ Cell order:
# ╠═0ba701f0-9006-11eb-08c9-6183e18fa170
# ╠═25a5f9d0-9006-11eb-3c3b-eb70551516f0
# ╠═91aa6940-9006-11eb-357e-c5e62eee94cd
# ╠═a7403aa0-9006-11eb-2825-3d0ac0de4494
# ╠═af2bd120-9006-11eb-2eb3-87c224af2b6a
# ╠═b5562fa0-9006-11eb-0e28-cdf834acd9e7
# ╠═5d4d630e-0a68-415b-a9cd-8650c6875914
# ╠═c711a7d4-af33-4445-823b-514c1dbcbf7b
# ╟─f5a8aced-e91e-4b63-9eeb-8f7d32d5dc9b
# ╟─076fdfdf-096d-4b66-ae81-081d5864c1ec
# ╟─5a3553a3-4af0-4ef3-b8b6-a7e635184636
