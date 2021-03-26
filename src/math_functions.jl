#using FFTW
using ImageFiltering

@. gaussian_model(x, (x₀, σ, A, B)) = A*exp(-(x - x₀)^2/(2σ^2)) + B
@. linear_model(x, (m, b)) = m*x + b

function convolve1d(yy, kernel)
    imfilter(yy, reflect(centered(kernel)), Fill(0))
    # NOTE: why? perhaps the Computational Thinking course will help...
end