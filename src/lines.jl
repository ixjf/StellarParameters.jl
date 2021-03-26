using Plots
using Peaks
using LsqFit
using MonteCarloMeasurements
using LinearAlgebra

function central_moving_average_filter(N, yy)
    for i=1:N
        yy = [yy[1], yy..., yy[end]]
        yy = [(y₁ + y₂ + y₃)/3 for (y₁, y₂, y₃) in zip(yy[begin:end-2], yy[begin+1:end-1], yy[begin+2:end])]
        # not the most optimized, but good enough
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
function isolate_line((λ, flux), λ_c_approx) # assumes line exists
    λ_interval = @. λ_c_approx - 0.5 <= λ <= λ_c_approx + 0.5 # Å
    λ = λ[λ_interval]
    flux = flux[λ_interval]

    flux_filtered = central_moving_average_filter(SPECTRUM_FILTER_ITERATIONS, flux)
    
    left_bound, right_bound = find_first_peaks_around(
        λ, 
        flux_filtered, 
        λ_c_approx, 
        SPECTRUM_PEAK_WINDOW_SIZE)

    λ = λ[left_bound:right_bound]
    flux = flux[left_bound:right_bound]

    (λ, flux)
end

function fit_line_as_gaussian((λ, flux), λ_c_corrected)
    x₀, σ, A, B = nothing, nothing, nothing, nothing
    x₀_error, σ_error, A_error, B_error = nothing, nothing, nothing, nothing

    #println(λ_c_corrected)

    try
        fit = curve_fit(gaussian_model, λ, flux, [λ_c_corrected, 0.1, -1.0, 1.0])

        !fit.converged && return nothing #begin; println("!converged"); return nothing; end;

        x₀, σ, A, B = coef(fit)

        x₀_error, σ_error, A_error, B_error = stderror(fit)
    catch e
        #scatter(λ, flux; legend=false)
        #gui()
        #sleep(0.2)
        #println("caught e: $e")
        #println(stacktrace(catch_backtrace()))
        return nothing
    end
 
    x₀ = x₀ ± x₀_error
    σ = σ ± σ_error
    A = A ± A_error
    B = B ± B_error

    W_λ = calculate_equivalent_width(σ, A, B)

    if mean(W_λ) < 0 || std(W_λ)/mean(W_λ)*100 > W_λ_UNC_THRESHOLD
        # scatter(λ, flux)
        #println(length(gaussian_model.(λ, (x₀, σ, A, B))))
        # scatter!(λ, gaussian_model(λ, (mean(x₀), mean(σ), mean(A), mean(B))); legend=false)
        # gui()
        #println("bad W_λ")
        return nothing
    end

    (x₀, σ, A, B)
end

function calculate_equivalent_width(σ, A, B)
    -sqrt(2π)*A*σ/B
end

function correct_λ_for_vrad(λ_c, vrad)
    λ_c*(1 + vrad/(c*10^(-3)))
end