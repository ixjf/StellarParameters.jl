# TODO: should this be a package? or an independent command line tool?

#using Profile

#import Measurements: ±, value, uncertainty

include("config.jl")

#±(μ::Real, σ) = μ + σ*Particles(MONTE_CARLO_NUM_SAMPLES)
#±(μ::Missing, σ::Missing) = missing + missing*Particles(MONTE_CARLO_NUM_SAMPLES)

include("thread_safe_dict.jl")
include("math_functions.jl")
include("spectra.jl")
include("line_list.jl")
include("excitation_temp.jl")
include("lines.jl")
include("profiles.jl")
include("grid_compare.jl")

# use @avx everywhere (LoopVectorization)
# parallelize (no methods tried worked well; multithreaded always turns out to be
# slower - isolate_all_lines_found_in_spectrum is slow af)
# optimize memory usage
# optimize monte carlo measurements (see performance tips)

# document code

# 1. config should be overrideable
# 2. all the filter()ing in excitation temp still makes me a bit crazy
# 3. other than that, finish writing main alg, make sure everything works
# then do analysis
# + optimization at same time
# 4. grid cache define type

# begin
#     full_line_list, ew_list_Sun = read_line_list()
#     multiplet_list = group_lines_into_multiplets(full_line_list)

#     #ew_cache = ThreadSafeDict{Tuple{String, Int64, Float64}, Dict{Float64, Particles{Float64, MONTE_CARLO_NUM_SAMPLES}}}()
#     grid_cache = GridCache()

#     sisma_archive_path = joinpath(@__DIR__, "../../data/sisma_archive/")
#     for filename in filter(x -> endswith(x, "fits"), readdir(sisma_archive_path))
#         λ, flux, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))

#         vrad = opt_data["VRAD"]

#         lines_obs, ew_list_obs = isolate_all_lines_found_in_spectrum(
#             vrad, 
#             nothing,
#             collect(keys(full_line_list)), 
#             (λ, flux))
        
#         best_combo_i, best_combo_j, _ = find_best_multiplet_combo_for_Texc_estimate(
#             multiplet_list,
#             lines_obs,
#             ew_list_Sun)
            
#         _, _, Texc = find_best_multiplet_combo_for_Texc_estimate(
#             [multiplet_list[best_combo_i], multiplet_list[best_combo_j]],
#             lines_obs,
#             ew_list_obs)
                
#         Teff = opt_data["TEFF"]
#         println("$filename: best combo: ($best_combo_i, $best_combo_j) Texc: $Texc Teff: $Teff")
    
#         R = opt_data["R"]
#         vsini = opt_data["VSINI"]
        
#         Teff, logg, FeH, α_Fe = match_ew_against_grid(ew_list_obs, Texc, R, vsini, grid_cache)#, ew_cache)
#         # println("match_ew_against_grid")

#         println("Multiplet: ($(best_combo_i), $(best_combo_j)) \t Texc: $Texc \t Teff: $Teff \t log(g): $logg \t [Fe/H]: $FeH \t [α/Fe]: $α_Fe")
#     end
# end

# full_line_list, ew_list_Sun = read_line_list()

# ###

# multiplet_list = group_lines_into_multiplets(full_line_list)

# ###

# sisma_archive_path = joinpath(@__DIR__, "../../data/sisma_archive/")

# ###

# d = Dict()

# #i1 = 0
# for filename in filter(x -> endswith(x, "fits"), readdir(sisma_archive_path))
#     #global i1 += 1

#     #i1 != 27 && continue
    
#     #filename != "HD169822_20120724_0000_nor.fits" && continue

#     λ, flux, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))
    
#     #opt_data["TEFF"] != 5700 && continue

#     vrad = opt_data["VRAD"]
#     lines, ew_list = isolate_all_lines_found_in_spectrum(
#         vrad,
#         nothing,
#         keys(full_line_list),
#         (λ, flux))
    
#     push!(d, filename => (λ, flux, opt_data, lines, ew_list))
# end
    
# ###

# p = scatter([], [])

# min_Teff = +Inf
# max_Teff = -Inf

# grid_cache = GridCache()

# for (filename, (λ, flux, opt_data, lines, ew_list)) in d
#     best_combo_i, best_combo_j, _ = find_best_multiplet_combo_for_Texc_estimate(
#         multiplet_list, 
#         lines, 
#         ew_list_Sun)

#     _, _, Texc = find_best_multiplet_combo_for_Texc_estimate(
#         [multiplet_list[best_combo_i], multiplet_list[best_combo_j]],
#         lines,
#         ew_list)
        
#     Teff = opt_data["TEFF"]
    
#     println("best combo: ($best_combo_i, $best_combo_j) Texc: $Texc Teff: $Teff")

#     scatter!([Teff], [Texc]; legend=false)
    
#     global min_Teff = min(Teff, min_Teff)
#     global max_Teff = max(Teff, max_Teff)

#     Teff, logg, FeH, α_Fe = match_ew_against_grid(ew_list, Texc, opt_data["R"], opt_data["VSINI"], grid_cache)

#     println("Teff: $Teff (ref: $(opt_data["TEFF"])) logg: $logg (ref: $(opt_data["LOGG"])), FeH: $FeH (ref: $(opt_data["FEH"])) α_Fe: $α_Fe")
# end

# plot!([min_Teff, max_Teff], [min_Teff, max_Teff]; legend=false)