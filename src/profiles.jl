using Statistics

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