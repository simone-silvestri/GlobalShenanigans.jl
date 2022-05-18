using JLD2

global_filepath = "/storage1/uq-global/GlobalShenanigans.jl/"
ic_filepath = global_filepath * "quick_check.jld2"
println("loading from ", ic_filepath)
jlfile = jldopen(ic_filepath)

tkeys = keys(jlfile["timeseries"]["T"])[3:end]
# 1:12 first year, 1+1*12:2*12 second year, 1+2*12:3*12 third year and so forth
Ttrac = jlfile["timeseries"]["T"][tkeys[1]]
Strac = jlfile["timeseries"]["S"][tkeys[1]]
Ttrac[continentst] .= NaN
Strac[continentst] .= NaN