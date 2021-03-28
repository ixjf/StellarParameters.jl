using Plots
using LsqFit
using Statistics
using Measurements: value, uncertainty, ±

function fit_curve_of_growth(xx, yy)
    p0 = [1.0, 0.0]
    local fit

    while true
        fit = curve_fit(linear_model, xx, yy, p0)

        # # TODO: is this necessary? looks like all fits are linear!
        # sum_residuals = sum(residuals(fit))
        # if sum_residuals > MULTIPLET_LINFIT_RESIDUALS_THRESHOLD
        #     return nothing
        # end

        σ = √(mse(fit))
        m, b = coef(fit)
        yy_pred = [m*x + b for x in xx]

        outliers = []

        for i=1:length(yy)
            y = yy[i]
            x = xx[i]
            y_pred = yy_pred[i]

            if abs(y - y_pred) > σ*CURVE_OF_GROWTH_SIGMACLIP
                push!(outliers, i)
            end
        end

        if length(outliers) == 0
            break
        end

        deleteat!(yy, outliers)
        deleteat!(xx, outliers)
    end

    fit
end

function calculate_curve_of_growth(multiplet, lines_used, ew_list)
    xx, yy = [], [] # TODO: make yy type stable

    for (λ_c_approx, line_data) in filter(kv -> kv[1] in keys(lines_used), multiplet)
        λ_c_exact = value(lines_used[λ_c_approx][1])
        push!(xx, line_data.loggf + log10(λ_c_exact))
        push!(yy, log10(value(ew_list[λ_c_approx])/λ_c_exact))
        # TODO: should I propagate uncertainty of EW to here?
        # this is a problem because some of the particles in the particle
        # distribution for some W_λ might have W_λ negative, so
        # log10 comes out wrong.
    end

    (xx, yy)
end

function group_lines_into_multiplets(line_list) # TODO: take multiplet EP ranges as input
    multiplets = []
    
    for bounds in MULTIPLETS
        multiplet = filter(kv -> bounds[1] <= kv[2].χ <= bounds[2], line_list)
        push!(multiplets, multiplet)
    end

    multiplets
end

function find_best_multiplet_combo_for_Texc_estimate(alg, multiplet_list, lines_used, ew_list)
    prev_error = +Inf # how exact (compared to the Sun's temperature) it is
    best_combo_i, best_combo_j = 0, 0
    Texc = 0.0 ± 0.0

    for i=1:length(multiplet_list), j=i+1:length(multiplet_list)
        xx1, yy1 = calculate_curve_of_growth(multiplet_list[i], lines_used, ew_list)
        fit1 = fit_curve_of_growth(xx1, yy1)

        # TODO: same as below
        multiplet_i = filter(kv -> kv[1] in keys(lines_used), multiplet_list[i])
        χ₁ = sum([line_data.χ for (_, line_data) in multiplet_i])/length(multiplet_i)

        xx2, yy2 = calculate_curve_of_growth(multiplet_list[j], lines_used, ew_list)
        fit2 = fit_curve_of_growth(xx2, yy2)

        # TODO: Doing this multiple times for the same multiplets...
        multiplet_j = filter(kv -> kv[1] in keys(lines_used), multiplet_list[j])
        χ₂ = sum([line_data.χ for (_, line_data) in multiplet_j])/length(multiplet_j)

        new_Texc_est = calculate_Texc_from_COG_pair(alg, (χ₁, xx1, fit1), (χ₂, xx2, fit2))

        let new_error = abs(value(new_Texc_est) - SUN_REF_TEFF)
            if (new_error < prev_error) || (new_error == prev_error && uncertainty(new_Texc_est) < uncertainty(Texc))
                best_combo_i, best_combo_j = i, j
                prev_error = new_error
                Texc = new_Texc_est
            end
        end

        # scatter(xx1, mean.(yy1); yerror=std.([ew_list[x] for (x,_) in multiplet_i]))
        # plot!(xx1, [fit1.param[1]*x + fit1.param[2] for x in xx1])
        # scatter!(xx2, mean.(yy2); yerror=std.([ew_list[x] for (x,_) in multiplet_j]))
        # plot!(xx2, [fit2.param[1]*x + fit2.param[2] for x in xx2])
        # gui()
        # sleep(2)
    end

    (best_combo_i, best_combo_j, Texc)
end

function weighted_median(m1, b1, xx1, m2, b2, xx2)
    y_values = vcat(
        [m1*x + b1 for x in xx1], 
        [m2*x + b2 for x in xx2])

    y_median = median(y_values) # median is least affected by outliers
    # so it's the best property to describe a weighted central tendency
end

function unweighted_median(m1, b1, xx1, m2, b2, xx2)
    y_min_multiplet1, y_max_multiplet1 = (m1*xx1[begin] + b1, m1*xx1[end] + b1)
    y_min_multiplet2, y_max_multiplet2 = (m2*xx2[begin] + b2, m2*xx2[end] + b2)

    y_min = max(y_min_multiplet1, y_min_multiplet2)
    y_max = min(y_max_multiplet1, y_max_multiplet2)

    # TODO: 1. weighted median is too affected by outliers.
    # outliers should be removed (Texc is wild!)

    y_median = (y_min + y_max)/2
end

# NOTE: unweighted_median: this method has the problem of there having to be points
# defined in the same interval of y for the two multiplets
# NOTE: make alg a separate function. so you pass this function as an argument.
function calculate_Texc_from_COG_pair(alg, (χ_multiplet1, xx_multiplet1, fit_multiplet1), 
        (χ_multiplet2, xx_multiplet2, fit_multiplet2))
    m_multiplet1_err, b_multiplet1_err = stderror(fit_multiplet1)
    m_multiplet2_err, b_multiplet2_err = stderror(fit_multiplet2)

    m_multiplet1, b_multiplet1 = coef(fit_multiplet1)
    m_multiplet1 = m_multiplet1 ± m_multiplet1_err
    b_multiplet1 = b_multiplet1 ± b_multiplet1_err

    m_multiplet2, b_multiplet2 = coef(fit_multiplet2)
    m_multiplet2 = m_multiplet2 ± m_multiplet2_err
    b_multiplet2 = b_multiplet2 ± b_multiplet2_err

    y_median = alg(m_multiplet1, b_multiplet1, xx_multiplet1, 
                   m_multiplet2, b_multiplet2, xx_multiplet2)

    x_multiplet1_median = (y_median - b_multiplet1)/m_multiplet1
    x_multiplet2_median = (y_median - b_multiplet2)/m_multiplet2

    Δ = abs(x_multiplet1_median - x_multiplet2_median)

    Δχ = abs(χ_multiplet1 - χ_multiplet2)

    Texc = 5040*Δχ/Δ
end