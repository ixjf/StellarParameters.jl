using DelimitedFiles

Symbol("FeI")
Symbol("FeII")

struct AtomicParameters
    El::Symbol
    χ::Float64
    loggf::Float64
end

function read_line_list()
    data = readdlm(LINE_LIST_FILEPATH; comments=true, comment_char='#', header=false)

    line_list = Dict{Float64, AtomicParameters}()
    ew_list = Dict{Float64, Float64}()

    for row in eachrow(data)
        λ_c_approx = row[1]
        χ = row[2]
        loggf = row[3]
        El = row[4] == "FeI" ? :FeI : :FeII
        W_λ_Sun = row[5]*10^(-3) # Å
        push!(line_list, λ_c_approx => AtomicParameters(El, χ, loggf))
        push!(ew_list, λ_c_approx => W_λ_Sun)
    end

    (line_list, ew_list)
end