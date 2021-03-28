#using MonteCarloMeasurements: mean, std, @unsafe
using Measurements: Measurement, value, uncertainty, ±
#using LsqFit
#using FFTW
#using Plots
#using Base.Threads
using Folds

# NOTE: move spectra data to a class?

GridCache = ThreadSafeDict{String, Array{Float64, 1}}

function convolve_simulated_profile(R, vsini, λ_c, continuum_flux, λ, flux)
    flux = flux[:] # copy, don't want to alter `flux`
    flux .-= value(continuum_flux) # remove constant factor before convolution

    if !ismissing(vsini)
        Δλ = λ[2] - λ[1]
        rp_kernel = rotational_profile(Δλ, EPSILON, value(λ_c), value(vsini))
        flux = convolve1d(flux, rp_kernel)
    end

    ip_kernel = instrumental_profile((λ, flux), R)
    flux = convolve1d(flux, ip_kernel)

    flux .+= value(continuum_flux)

    (λ, flux)
end

function isolate_all_lines_found_in_spectrum(vrad, profile_fn, line_list, (λ, flux)) # line_list = tsantaki's
    lines = Dict{Float64, Tuple{Measurement{Float64}, Array{Float64}, Array{Float64}}}()
    #ew_list = Dict{Float64, Particles}()
    ew_list = Dict{Float64, Measurement{Float64}}()
    # TODO: make previous types type-safe

    for λ_c_approx in line_list
        !(λ[begin] <= λ_c_approx <= λ[end]) && continue

        λ_c_corrected = !ismissing(vrad) ? correct_λ_for_vrad(λ_c_approx, value(vrad)) : λ_c_approx

        λ_line, flux_line = isolate_line((λ, flux), λ_c_corrected)

        fit_result = fit_line_as_gaussian((λ_line, flux_line), λ_c_corrected)

        isnothing(fit_result) && continue

        λ_c_exact, σ, A, B = fit_result
        
        if profile_fn !== nothing
            # B = continuum_flux
            #scatter(λ_line, flux_line)
            λ_line, flux_line = profile_fn(λ_c_corrected, B, λ_line, flux_line)
            #scatter!(λ_line, flux_line; legend=false)
            #gui()
            #sleep(0.2)
            #println(mean(λ_c_exact))

            fit_result = fit_line_as_gaussian((λ_line, flux_line), value(λ_c_exact))

            isnothing(fit_result) && continue
            # NOTE: the above may fail (convolution broadens the line too much
            # and LsqFit can't find a good fit?)

            λ_c_exact, σ, A, B = fit_result
        end

        W_λ = calculate_equivalent_width(σ, A, B)

        push!(lines, λ_c_approx => (λ_c_exact, λ_line, flux_line))
        push!(ew_list, λ_c_approx => W_λ)
    end

    (lines, ew_list)
end

function calculate_min_squared_error_ew(ew_list_obs, ew_list_grid)
    found_in_both = intersect(keys(ew_list_obs), keys(ew_list_grid))
    lines_to_compare = filter(x -> x[1] in found_in_both, ew_list_obs)

    #χ² = 0.0 + 0.0*Particles(MONTE_CARLO_NUM_SAMPLES)
    χ² = 0.0 ± 0.0

    if length(lines_to_compare) == 0
        #return value(χ²) ± uncertainty(χ²)
        return χ²
    end

    for λ_c in keys(lines_to_compare)
        W_λ_match = ew_list_obs[λ_c]
        W_λ_grid = ew_list_grid[λ_c]

        #W_λ_match = measurement(mean(ew_list_obs[λ_c]), std(ew_list_obs[λ_c]))
        #W_λ_grid = measurement(mean(ew_list_grid[λ_c]), std(ew_list_grid[λ_c]))

        χ² += (W_λ_match - W_λ_grid)^2
    end
    
    χ² /= length(lines_to_compare)
    
    #println(χ²)

    χ²
    #value(χ²) ± uncertainty(χ²)
end

# TODO: can I reuse EWs from the grid? If same R and vsini.
function match_ew_against_grid(ew_list_obs, Texc, R, vsini, grid_cache)#, ew_cache)
    λ_grid = get!(grid_cache, GRID_SPECTRA_WAVE_PATH) do
        read_spectrum_ambre_grid(GRID_SPECTRA_WAVE_PATH)
    end

    spectra = Dict{Tuple{Int64, Float64, Float64, Float64}, Array{Float64, 1}}()
    
    for filename in readdir(GRID_SPECTRA_FLUX_PATH)
        Teff = parse(Int64, filename[2:5])
        
        !(value(Texc) - uncertainty(Texc) <= Teff <= value(Texc) + uncertainty(Texc)) && continue
        
        logg = parse(Float64, filename[8:11])
        FeH = parse(Float64, filename[23:27])
        α_Fe = parse(Float64, filename[30:34])
        
        filepath = joinpath(GRID_SPECTRA_FLUX_PATH, filename)
        
        flux_grid = get!(grid_cache, filepath) do
            read_spectrum_ambre_grid(filepath)
        end
        
        push!(spectra, (Teff, logg, FeH, α_Fe) => flux_grid)
    end
    
    # matches = ThreadSafeDict{Tuple{Int64, Float64, Float64, Float64}, Particles{Float64, MONTE_CARLO_NUM_SAMPLES}}()
    
    matches = Folds.mapreduce(merge, spectra; init=Dict()) do (parameters_grid, flux_grid)
        #println("\t Matching against: $filename")

        #ew_list_grid = get!(ew_cache, (filename, R, !ismissing(vsini) ? mean(vsini) : 0.0)) do
            # println("isolate_all_lines_found_in_spectrum")
            lines_grid, ew_list_grid = isolate_all_lines_found_in_spectrum(
                0.0, 
                (λ_c, continuum_flux, λ_line, flux_line) -> convolve_simulated_profile(R, vsini, λ_c, continuum_flux, λ_line, flux_line),
                keys(ew_list_obs), 
                (λ_grid, flux_grid))
            # FIXME: the above is sometimes taking >9 seconds!
        #    ew_list_grid
        #end

        χ² = calculate_min_squared_error_ew(ew_list_obs, ew_list_grid)

        Dict(parameters_grid => χ²)
        # push!(matches, parameters_grid => χ²)
    end

    χ²_final, (Teff_final, logg_final, FeH_final, α_Fe_final) = findmin(matches)#@unsafe findmin(matches)

    ΔTeff, Δlogg, ΔFeH, Δα_Fe = 0.0, 0.0, 0.0, 0.0

    for ((Teff, logg, FeH, α_Fe), χ²) in matches
        # does the uncertainty interval for χ² intersect the uncertainty
        # interval for χ²? 
        if value(χ²_final) - uncertainty(χ²_final) <= value(χ²) + uncertainty(χ²) &&
           value(χ²_final) + uncertainty(χ²_final) >= value(χ²) - uncertainty(χ²)
            ΔTeff = max(abs(Teff_final - Teff), ΔTeff)
            Δlogg = max(abs(logg_final - logg), Δlogg)
            ΔFeH = max(abs(FeH_final - FeH), ΔFeH)
            Δα_Fe = max(abs(α_Fe_final - α_Fe), Δα_Fe)
        end
    end

    (Teff_final ± ΔTeff, 
     logg_final ± Δlogg, 
     FeH_final ± ΔFeH, 
     α_Fe_final ± Δα_Fe)
end