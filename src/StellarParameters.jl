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
# 7. (TODO) Preallocate memory so we don't get absurd memory usage.
# Also since no profile stuff in match_ew_against_grid anymore, I can easily cache 
# the results of comparison against one spectrum.

# document code

# 1. config should be overrideable
# 2. all the filter()ing in excitation temp still makes me a bit crazy
# + optimization at same time

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

# # ###

# sisma_archive_path = joinpath(@__DIR__, "../../data/sisma_archive/")

# # ###

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
#         keys(full_line_list),
#         (λ, flux))
    
#     push!(d, filename => (λ, flux, opt_data, lines, ew_list))
# end
    
# ###

# p = scatter([], [])

# min_Teff = +Inf
# max_Teff = -Inf

# for filename in filter(x -> endswith(x, "fits"), readdir(sisma_archive_path))	
# 	λ, flux, opt_data = read_spectrum_sisma(joinpath(sisma_archive_path, filename))
	
# 	vrad = opt_data["VRAD"]
	
# 	lines, ew_list = isolate_all_lines_found_in_spectrum(
# 		vrad,
# 		keys(full_line_list),
# 		(λ, flux),
# 		5,
# 		5,
# 		12.5)
	
# 	grid_cache = GridCache()
	
# 	#best_combo_i, best_combo_j, _ = find_best_multiplet_combo_for_Texc_estimate(
# 	#	weighted_median,
# 	#	multiplet_list, 
# 	#	Dict([λ_c_approx => (λ_c_approx,) for (λ_c_approx,_) in lines]), 
# 	#	ew_list_Sun)

# 	_, _, Texc = find_best_multiplet_combo_for_Texc_estimate(
# 		weighted_median,
# 		[multiplet_list[1], multiplet_list[4]],
# 		lines,
# 		ew_list)

# 	println("--------------------------")
	
# 	#Teff, logg, FeH, α_Fe = match_ew_against_grid(ew_list, Texc, opt_data["R"], opt_data["VSINI"], grid_cache)

#     #println("Teff: $Teff logg: $logg FeH: $FeH α_Fe: $α_Fe")
# end

# grid_cache = GridCache()

# for (filename, (λ, flux, opt_data, lines, ew_list)) in d
#     best_combo_i, best_combo_j, _ = find_best_multiplet_combo_for_Texc_estimate(
#         weighted_median,
#         multiplet_list, 
#         Dict([λ_c_approx => (λ_c_approx,) for (λ_c_approx,_) in lines]), 
#         ew_list_Sun)

#     _, _, Texc = find_best_multiplet_combo_for_Texc_estimate(
#         weighted_median,
#         [multiplet_list[best_combo_i], multiplet_list[best_combo_j]],
#         lines,
#         ew_list)
    
#     Teff = opt_data["TEFF"]
    
#     #println("best combo: ($best_combo_i, $best_combo_j) Texc: $Texc Teff: $Teff")

#     #scatter!([Teff], [Texc]; legend=false)
    
#     #global min_Teff = min(Teff, min_Teff)
#     #global max_Teff = max(Teff, max_Teff)

#     Teff, logg, FeH, α_Fe = match_ew_against_grid(ew_list, Texc, opt_data["R"], opt_data["VSINI"], grid_cache)

#     println("Teff: $Teff (ref: $(opt_data["TEFF"])) logg: $logg (ref: $(opt_data["LOGG"])), FeH: $FeH (ref: $(opt_data["FEH"])) α_Fe: $α_Fe")
# end

# plot!([min_Teff, max_Teff], [min_Teff, max_Teff]; legend=false)

# begin
# 	results = []
	
# 	grid_cache = GridCache()
	
# 	for filename in ["estrela2_vcor.fits"]
# 		λ, flux, Δx = read_spectrum_unknown_source_prof(joinpath(@__DIR__, "../../data/estrelas_a_analisar/", filename))
		 
# 		lines, ew_list = isolate_all_lines_found_in_spectrum(
# 			0.0,
# 			keys(full_line_list),
# 			(λ, flux),
# 			5, # already determined in 2_lines.jl
# 			5, # already determined in 2_lines.jl
# 			12.5) # already determined in 2_lines.jl

# 		best_combo_i, best_combo_j, _ = find_best_multiplet_combo_for_Texc_estimate(
# 			weighted_median, 
# 			multiplet_list, 
# 			Dict([λ_c_approx => (λ_c_approx,) for (λ_c_approx,_) in lines]), 
# 			ew_list_Sun,
# 			2)

# 		_, _, Texc = find_best_multiplet_combo_for_Texc_estimate(
# 			weighted_median, 
# 			[multiplet_list[best_combo_i], multiplet_list[best_combo_j]],
# 			lines,
# 			ew_list,
# 			2)

# 		params = match_ew_against_grid(ew_list, Texc, grid_cache)
		
# 		isnothing(params) && continue

# 		Teff, logg, FeH, α_Fe = params

# 		push!(
# 			results, 
# 			filename => (λ, 
# 				flux, 
# 				lines, 
# 				ew_list, 
# 		      	best_combo_i,
# 		      	best_combo_j,
# 				Texc, 
# 				Teff, 
# 				logg, 
# 				FeH, 
# 				α_Fe))
# 	end
# end

# f = FITS("../data/estrelas_a_analisar/estrela2_vcor.fits")
# hdu = f[1]

# header = read_header(hdu)

# println(header)

# TODO: this is a mess! multiplet_list here, lines_used, ew_list! I don't
# know which λ_c the function will use, and it's not obvious that lines
# must have λ_c_exact = λ_c_approx for the Sun but that the λ and flux are not required
# in the tuple!