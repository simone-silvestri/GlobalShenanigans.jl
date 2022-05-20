using JLD2

function load_statistics(; ic_filepath = "/storage1/uq-global/GlobalShenanigans.jl/")

    println("loading from ", ic_filepath)
    jlfile = jldopen(ic_filepath)

    tkeys = keys(jlfile["timeseries"]["T"])[3:end]

    # 1:12 first year, 1+1*12:2*12 second year, 1+2*12:3*12 third year and so forth
    Ttrac = zeros(size(jlfile["timeseries"]["T"][tkeys[end]]))
    Strac = zeros(size(jlfile["timeseries"]["S"][tkeys[end]]))

    for key in tkeys[1:24]
        Ttrac .+= jlfile["timeseries"]["T"][key] / 24
        Strac .+= jlfile["timeseries"]["S"][key] / 24
    end
    
    continentst = (Ttrac .== 0.0)
    continentss = (Strac .== 0.0)

    Ttrac[continentst] .= NaN
    Strac[continentss] .= NaN

    return Ttrac, Strac, tkeys
end

function load_single_file(ic_filepath = "/storage1/uq-global/GlobalShenanigans.jl/")
    println("loading from ", ic_filepath)
    jlfile = jldopen(ic_filepath)


    # 1:12 first year, 1+1*12:2*12 second year, 1+2*12:3*12 third year and so forth
    Ttrac = jlfile["u/data"][5:end-4, 5:end-4, 5:end-4]
    Strac = jlfile["v/data"][5:end-4, 5:end-4, 5:end-4]
    
    continentst = (Ttrac .== 0.0)
    continentss = (Strac .== 0.0)

    Ttrac[continentst] .= NaN
    Strac[continentss] .= NaN

    return Ttrac, Strac
end
