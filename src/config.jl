const c = 3e8 # m s-1

const LINE_LIST_FILEPATH = joinpath(@__DIR__, "../../data/line_list_tsantaki.dat")
const MULTIPLETS = [(2.18, 2.28), (2.42, 2.61), (3.57, 3.69), (4.19, 4.22), (4.55, 4.65)]
# const MULTIPLET_LINFIT_RESIDUALS_THRESHOLD = 1e-7
#const EXCITATION_TEMP_ALG = :weighted_median
const SPECTRUM_FILTER_ITERATIONS = 3
const SPECTRUM_PEAK_WINDOW_SIZE = 5
#const SPECTRUM_FILTER_CONTINUUM_ITERATIONS = 10000
const GRID_SPECTRA_FLUX_PATH = joinpath(@__DIR__, "../../data/ambre_proj/GES_UVESRed580_deltaAlphaFe+0.0_fits/")
const GRID_SPECTRA_WAVE_PATH = joinpath(@__DIR__, "../../data/ambre_proj/GES_UVESRed580_Lambda.fits")
const EPSILON = 0.6
const SUN_REF_TEFF = 5_777
const W_λ_UNC_THRESHOLD = 10 # %
#const MONTE_CARLO_NUM_SAMPLES = 30_000
const CURVE_OF_GROWTH_SIGMACLIP = 2



# 1. which multiplets were used?
# 2. is unweighted_median better?