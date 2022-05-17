using NCDatasets, Statistics, HDF5

include("ecco_wrangle.jl")

# File Format Stuff
field, latitude, longitude, raditude = grab_ecco_field("THETA", 1, 1)
ecco_size = size(field)
# oceananigans_format = zeros((ecco_size[2], ecco_size[1], ecco_size[3], 12))
oceananigans_format = zeros((ecco_size[1], ecco_size[2], ecco_size[3], 12))

fname = "wrangle_ecco_3.h5"
fid = h5open(fname, "w")
create_group(fid, "fields")
create_group(fid, "grid")
fid["grid"]["latitude"] = latitude
fid["grid"]["longitude"] = longitude
fid["grid"]["raditude"] = raditude

# 1:24, ..., 3 * 4 * 2, we have 26 years but this is good enough
year_groups = [4*(i-1)+1:4*(i-1)+4 for i in 1:6]
year_groups = [1:5] # for now only handle one average

month_integer = 1
month_dictionary = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
fieldnames = ["THETA", "SALT"]

for fieldname in fieldnames
    println("Looking at field ", fieldname)
    for month_integer in 1:12
        println("Aggregating month ", month_dictionary[month_integer])

        averaged_fields = []
        for years in year_groups
            println("Currently on year group ", years)
            fields = []
            for year_integer in years
                println("At year ", year_integer + 1991)
                field, _, _, _ = grab_ecco_field(fieldname, year_integer, month_integer)
                push!(fields, field)
                quantile_lims = (0.01, 0.99)
                quantile_field = round.(Ref(Int), quantile.(Ref((field[(!).(isnan.(field))])), quantile_lims))
                println("the rounded ", quantile_lims, " quantiles of ", fieldname, " are ", quantile_field)
            end
            println("averaging")
            push!(averaged_fields, mean(fields))
            println("--------------------------------")
        end

        ecco_format = averaged_fields[1]
        # change_format!(view(oceananigans_format, :, :, :, month_integer), ecco_format)
        view_oceananigans_format = view(oceananigans_format, :, :, :, month_integer)
        view_oceananigans_format .= ecco_format
    end

    if fieldname == "THETA"
        fid["fields"]["temperature"] = oceananigans_format
    elseif fieldname == "SALT"
        fid["fields"]["salinity"] = oceananigans_format
    end
end

close(fid)
