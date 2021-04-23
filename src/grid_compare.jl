#using MonteCarloMeasurements: mean, std, @unsafe
using Measurements: Measurement, value, uncertainty, ±
#using LsqFit
#using FFTW
#using Plots
#using Base.Threads
using Folds

# NOTE: move spectra data to a class?

# not multi threaded, unnecessary right now (for it to be thread-safe)
GridCache = ThreadSafeDict{String, Array{Float64, 1}}

function calculate_min_squared_error_ew(ew_list_obs, ew_list_grid)
    found_in_both = intersect(keys(ew_list_obs), keys(ew_list_grid))
    lines_to_compare = filter(x -> x[1] in found_in_both, ew_list_obs)

    #println("Comparing $(length(lines_to_compare)) lines")

    #χ² = 0.0 + 0.0*Particles(MONTE_CARLO_NUM_SAMPLES)
    χ² = 0.0 ± 0.0

    for λ_c in keys(lines_to_compare)
        W_λ_match = ew_list_obs[λ_c]
        W_λ_grid = ew_list_grid[λ_c]

        #W_λ_match = measurement(mean(ew_list_obs[λ_c]), std(ew_list_obs[λ_c]))
        #W_λ_grid = measurement(mean(ew_list_grid[λ_c]), std(ew_list_grid[λ_c]))

        χ² += (W_λ_match - W_λ_grid)^2

        #println(uncertainty((W_λ_match - W_λ_grid)^2))
    end
    
    χ² /= length(lines_to_compare)
    
    #println("Total χ²: $χ²")
    #println(χ²)

    χ²
    #value(χ²) ± uncertainty(χ²)
end

# TODO: can I reuse EWs from the grid? If same R and vsini.
function match_ew_against_grid(ew_list_obs, Texc, grid_cache)#, ew_cache)
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
    
    matches = Dict()
    
    for (parameters_grid, flux_grid) in spectra
        #println("\t Matching against: $filename")

        #ew_list_grid = get!(ew_cache, (filename, R, !ismissing(vsini) ? mean(vsini) : 0.0)) do
            # println("isolate_all_lines_found_in_spectrum")
            lines_grid, ew_list_grid = isolate_all_lines_found_in_spectrum(
                0.0, 
                keys(ew_list_obs), 
                (λ_grid, flux_grid),
                5, # TODO: make this param adjustable
                5, # TODO: same
                12.5) # TODO: same
            # FIXME: the above is sometimes taking >9 seconds!
        #    ew_list_grid
        #end

        χ² = calculate_min_squared_error_ew(ew_list_obs, ew_list_grid)

        push!(matches, parameters_grid => χ²)
    end

    isempty(matches) && return nothing

    χ²_final, (Teff_final, logg_final, FeH_final, α_Fe_final) = findmin(matches)#@unsafe findmin(matches)

    ΔTeff, Δlogg, ΔFeH, Δα_Fe = 0.0, 0.0, 0.0, 0.0

    for ((Teff, logg, FeH, α_Fe), χ²) in matches
        #println(χ²)

        # does the uncertainty interval for χ² intersect the uncertainty
        # interval for χ²? 
        if value(χ²_final) - uncertainty(χ²_final) <= value(χ²) + uncertainty(χ²) &&
           value(χ²_final) + uncertainty(χ²_final) >= value(χ²) - uncertainty(χ²)
            ΔTeff = max(abs(value(Teff_final - Teff)), ΔTeff)
            Δlogg = max(abs(value(logg_final - logg)), Δlogg)
            ΔFeH = max(abs(value(FeH_final - FeH)), ΔFeH)
            Δα_Fe = max(abs(value(α_Fe_final - α_Fe)), Δα_Fe)
        end
    end

    ΔTeff += GRID_TEFF_UNC
    Δlogg += GRID_LOGG_UNC
    ΔFeH += GRID_FEH_UNC
    Δα_Fe += GRID_α_FE_UNC

    (Teff_final ± ΔTeff, 
     logg_final ± Δlogg, 
     FeH_final ± ΔFeH, 
     α_Fe_final ± Δα_Fe)
end