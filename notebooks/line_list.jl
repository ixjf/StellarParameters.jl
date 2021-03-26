### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 576fff7e-8406-11eb-14bd-65f753dcd34e
using CSV, DataFrames

# ╔═╡ 68f83080-8409-11eb-3b14-d394416309b2
using Plots

# ╔═╡ fe9817b0-8407-11eb-1e68-cd353e5e0199
begin
	file = CSV.File("../../data/line_list_tsantaki.dat"; comment="#", header=false)
end

# ╔═╡ a7d76c90-8408-11eb-0bd2-65fa35584f48
FeI_lines = filter(r -> r[4] != "FeII", file)

# ╔═╡ 02c69ed0-840b-11eb-028c-0b9f41bd7a68
multiplets = [] # criteria are: isolated blocks with EP more or less constant (i.e. many points at pretty much the same EP, graph doesn't look like it's increasing). not really very precise, just eyeballing it.

# ╔═╡ 016943a0-8409-11eb-1b0c-a330a9e0098e
sort!(FeI_lines; by=r -> r[2])

# ╔═╡ 6c71c6e2-8409-11eb-3842-758b7b4e9cee
scatter(1:length(FeI_lines), [r[2] for r in FeI_lines]; legend=false)

# ╔═╡ 1a3125a2-840a-11eb-0463-af5177375be8
FeI_lines_2_3 = filter(r -> 2 <= r[2] <= 3, FeI_lines)

# ╔═╡ 8e76c59e-840a-11eb-1a5b-8bb5a6619b38
scatter(1:length(FeI_lines_2_3), [r[2] for r in FeI_lines_2_3]; legend=false)

# ╔═╡ 0c616bf0-840b-11eb-0f73-d54772f4bbda
push!(multiplets, filter(r -> 2.18 <= r[2] <= 2.28, FeI_lines))

# ╔═╡ 3b612120-840b-11eb-05a8-bd36838b6d64
push!(multiplets, filter(r -> 2.42 <= r[2] <= 2.61, FeI_lines))

# ╔═╡ 9483bb9e-840b-11eb-3cf0-71c9994770d2
FeI_lines_3_4 = filter(r -> 3 <= r[2] <= 4, FeI_lines)

# ╔═╡ b5b33760-840b-11eb-25b8-5927797a8f4c
scatter(1:length(FeI_lines_3_4), [r[2] for r in FeI_lines_3_4]; legend=false)

# ╔═╡ 0f40c8ae-840c-11eb-2634-bd56b5a249c8
push!(multiplets, filter(r -> 3.57 <= r[2] <= 3.69, FeI_lines))

# ╔═╡ 47672c6e-840c-11eb-1435-59a3ad276fec
FeI_lines_3_5 = filter(r -> 3.8 <= r[2] <= 5.1, FeI_lines)

# ╔═╡ 6668a4a0-840c-11eb-2f46-639c88ddc090
scatter(1:length(FeI_lines_3_5), [r[2] for r in FeI_lines_3_5]; legend=false)

# ╔═╡ 9bb95320-840c-11eb-02b7-6d2a54c8da9e
push!(multiplets, filter(r -> 4.19 <= r[2] <= 4.22, FeI_lines))

# ╔═╡ 86d32d40-840d-11eb-0957-77a72c913266
push!(multiplets, filter(r -> 4.55 <= r[2] <= 4.65, FeI_lines))

# ╔═╡ 699d1450-840f-11eb-14f0-31fbff75e950
# Writing multiplet list to a CSV file
begin
	const column_names = ["i", "λ / Å", "EP / eV", "loggf", "EW (Sun) / mÅ"]
	
	table = DataFrame(i=Int64[], λ=Float64[], EP=Float64[], loggf=Float64[], EW=Float64[])
	
	for (i, multiplet) in enumerate(multiplets)
		for row in multiplet
			push!(table, [i  row[1]  row[2]  row[3]  row[5]])
		end
	end
	
	CSV.write("../../data/line_list_analysed.csv", table; header=column_names, append=false)
end

# ╔═╡ Cell order:
# ╠═576fff7e-8406-11eb-14bd-65f753dcd34e
# ╠═fe9817b0-8407-11eb-1e68-cd353e5e0199
# ╠═a7d76c90-8408-11eb-0bd2-65fa35584f48
# ╠═02c69ed0-840b-11eb-028c-0b9f41bd7a68
# ╠═016943a0-8409-11eb-1b0c-a330a9e0098e
# ╠═68f83080-8409-11eb-3b14-d394416309b2
# ╠═6c71c6e2-8409-11eb-3842-758b7b4e9cee
# ╠═1a3125a2-840a-11eb-0463-af5177375be8
# ╠═8e76c59e-840a-11eb-1a5b-8bb5a6619b38
# ╠═0c616bf0-840b-11eb-0f73-d54772f4bbda
# ╠═3b612120-840b-11eb-05a8-bd36838b6d64
# ╠═9483bb9e-840b-11eb-3cf0-71c9994770d2
# ╠═b5b33760-840b-11eb-25b8-5927797a8f4c
# ╠═0f40c8ae-840c-11eb-2634-bd56b5a249c8
# ╠═47672c6e-840c-11eb-1435-59a3ad276fec
# ╠═6668a4a0-840c-11eb-2f46-639c88ddc090
# ╠═9bb95320-840c-11eb-02b7-6d2a54c8da9e
# ╠═86d32d40-840d-11eb-0957-77a72c913266
# ╠═699d1450-840f-11eb-14f0-31fbff75e950
