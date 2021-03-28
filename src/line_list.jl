using DelimitedFiles

struct AtomicParameters
    χ::Float64
    loggf::Float64
end

# TODO: use NamedTuples?
function read_line_list()
    data = readdlm(LINE_LIST_FILEPATH; comments=true, comment_char='#', header=false)

    data = data[data[:,4] .!= "FeII",:]

    line_list = Dict{Float64, AtomicParameters}()
    ew_list = Dict{Float64, Float64}()

    for row in eachrow(data[data[:,4] .!= "FeII",:])
        λ_c_approx = row[1]
        χ = row[2]
        loggf = row[3]
        W_λ_Sun = row[5]*10^(-3) # Å
        push!(line_list, λ_c_approx => AtomicParameters(χ, loggf))
        push!(ew_list, λ_c_approx => W_λ_Sun)
    end

    (line_list, ew_list)
end