using Statistics
using Measurements: value

function instrumental_profile((λ, flux), R)
    λ_mean = mean([λ[begin], λ[end]])
    FWHM = λ_mean/R
    σ = FWHM/(2*√(2*log(2)))
    Δx = λ[2] - λ[1]

    kernel_λ_space = -4*FWHM:Δx:4*FWHM
    kernel = gaussian_model(kernel_λ_space, (0.0, σ, 1.0, 0.0))
    kernel ./= sum(kernel) # normalize to [0,1]
    
    kernel
end

function rotational_profile(Δx, ϵ, λ_c, vsini) # vsini in km/s, Δx for spectrum to be
    # convolved
    vsini == zero(vsini) && throw(DomainError("vsini", "value cannot be zero"))
    
    Δλₘ = λ_c*(vsini*1e3)/c

    kernel_λ_space = -Δλₘ-1:Δx:Δλₘ+1
    
    term = max.((@. 1 - (kernel_λ_space/Δλₘ)^2), 0.0)

    kernel = @. (2*(1 - ϵ)*term^(1/2) + (π*ϵ/2)*term)/(π*Δλₘ*(1 - ϵ/3))
    kernel ./= sum(kernel) # normalize to [0,1]

    kernel
end

function convolve_simulated_profile(R, vsini, λ_c, continuum_flux, λ, flux)
    flux = flux[:] # copy, don't want to alter `flux`
    flux .-= value(continuum_flux) # remove constant factor before convolution

    # FIXME: convolve1d uses mode='same'. EW after convolution
    # shouldn't change, but it does because the line is truncated
    # and I can no longer measure the entire area.
    # This also means that I DON'T need to convolve with the simulated
    # profile to determine EW...
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