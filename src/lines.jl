using Plots
using Peaks
using LsqFit
using Measurements: ±, value, uncertainty, Measurement
using LinearAlgebra
using Statistics

function central_moving_average_filter(N, yy)
    for i=1:N
        yy = @inbounds [yy[1], yy..., yy[end]]
        yy = @inbounds [sum(@view yy[i:(i + 2)])/3 for i in 1:(length(yy) - 2)]
    end

    yy
end

function find_first_peaks_around(xx, yy, x_center, window_size)
    i_center = findmin(@. abs(xx - x_center))[2]

    peaks = argmaxima(yy, window_size)

    if length(peaks) == 0
        return (firstindex(xx), lastindex(xx))
    end

    sort!(peaks; by=x -> abs(x - i_center))

    if length(peaks) == 1
        if peaks[1] < i_center
            return (peaks[1], lastindex(xx))
        else
            return (firstindex(xx), peaks[1])
        end
    end

    left_peak = minimum(peaks[1:2])
    right_peak = maximum(peaks[1:2])

    (left_peak, right_peak)
end

# isolate_line creates a copy (slices are by-copy by default!)
function isolate_line((λ, flux), λ_c_approx, line_filter_iterations, line_filter_window_size) # assumes line exists
    λ_interval = @. λ_c_approx - 0.5 <= λ <= λ_c_approx + 0.5 # Å
    λ = λ[λ_interval]
    flux = flux[λ_interval]

    flux_filtered = central_moving_average_filter(line_filter_iterations, flux)
    
    left_bound, right_bound = find_first_peaks_around(
        λ, 
        flux_filtered, 
        λ_c_approx, 
        line_filter_window_size)

    λ = λ[left_bound:right_bound]
    flux = flux[left_bound:right_bound]

    (λ, flux)
end

function fit_line_as_gaussian((λ, flux), λ_c, W_λ_unc_threshold=W_λ_UNC_THRESHOLD)
    #println(λ_c_corrected)
    
    # TODO: improve initial guess => increase number of lines used
    p0 = [λ_c, 0.1, -1.0, mean(flux)]
    
    x₀, σ, A, B = nothing, nothing, nothing, nothing
    x₀_error, σ_error, A_error, B_error = nothing, nothing, nothing, nothing

    try
        # TODO: pass Jacobian
        fit = curve_fit(gaussian_model, λ, flux, p0)

        !fit.converged && return nothing #begin; println("!converged"); return nothing; end;
    
        x₀, σ, A, B = coef(fit)
    
        x₀_error, σ_error, A_error, B_error = stderror(fit)
    catch e
        # println("exception rip")
        #scatter(λ, flux; legend=false)
        #gui()
        #sleep(0.2)
        #println("caught e: $e")
        #println(stacktrace(catch_backtrace()))
        return nothing
    end

    # plot(λ, gaussian_model(λ, (x₀, σ, A, B)))
    # gui()
    # sleep(3)
 
    x₀ = x₀ ± x₀_error
    σ = σ ± σ_error
    A = A ± A_error
    B = B ± B_error

    W_λ = calculate_equivalent_width(σ, A, B)

    if W_λ < 0 || uncertainty(W_λ)/value(W_λ)*100 > W_λ_unc_threshold
        # scatter(λ, flux)
        #println(length(gaussian_model.(λ, (x₀, σ, A, B))))
        # scatter!(λ, gaussian_model(λ, (mean(x₀), mean(σ), mean(A), mean(B))); legend=false)
        # gui()
        #println("bad W_λ")
        return nothing
    end

    (x₀, σ, A, B, W_λ)
end

function calculate_equivalent_width(σ, A, B)
    W_λ = -sqrt(2π)*A*σ/B
    #value(W_λ) + uncertainty(W_λ)*Particles(MONTE_CARLO_NUM_SAMPLES)
end

function correct_λ_for_vrad(λ_c, vrad)
    λ_c*(1 + vrad/(c*10^(-3)))
end

function isolate_all_lines_found_in_spectrum(
        vrad, line_list, (λ, flux), # line_list = tsantaki's
        line_filter_iterations=SPECTRUM_FILTER_ITERATIONS, 
        line_filter_window_size=SPECTRUM_PEAK_WINDOW_SIZE,
        W_λ_unc_threshold=W_λ_UNC_THRESHOLD)
    lines = Dict{Float64, Tuple{Measurement{Float64}, Array{Float64}, Array{Float64}}}()
    #ew_list = Dict{Float64, Particles}()
    ew_list = Dict{Float64, Measurement{Float64}}()
    # TODO: make previous types type-safe

    for λ_c_approx in line_list
        λ_c_corrected = !ismissing(vrad) ? correct_λ_for_vrad(λ_c_approx, value(vrad)) : λ_c_approx
        
        if !(λ[begin] <= λ_c_corrected <= λ[end])
            #println("line not in spectrum [λ_c_corrected: $λ_c_corrected, λ_c_begin: $(λ[begin]), λ_c_end: $(λ[end])")
            continue
        end

        λ_line, flux_line = isolate_line(
            (λ, flux), 
            λ_c_corrected, 
            line_filter_iterations, 
            line_filter_window_size)
            
        fit_result = fit_line_as_gaussian((λ_line, flux_line), λ_c_corrected, W_λ_unc_threshold)
        
        #scatter(λ_line, flux_line)
        #gui()
        #sleep(1)

        if isnothing(fit_result)
            #println("nothing")
            continue
        end

        λ_c_exact, σ, A, B, W_λ = fit_result

        push!(lines, λ_c_approx => (λ_c_exact, λ_line, flux_line))
        push!(ew_list, λ_c_approx => W_λ)
    end

    #println(length(line_list), " ", i, " ", length(lines), " ", length(ew_list))

    (lines, ew_list)
end