using NCDatasets

month_int = 1
year_int = 1
fieldname = "SALT"
local_directory = "/mnt/podaac_drive/Version4/Release4/interp_monthly/"
skip = 1

# Find available years
years = readdir(local_directory * fieldname)
year = years[year_int]
# Find available months
println("Looking at year ", year)
months = readdir(joinpath(local_directory * fieldname, year))
month = months[month_int]
tic = time()
println("Looking at month ", month)
full_path = joinpath(local_directory, fieldname, year, month)
# Load Data Set
ds = Dataset(full_path, "r")
toc = time()
println("Time to load dataset = ", toc - tic, " seconds")

field = ds[fieldname]

# Plot variables
latitude = ds["latitude"][:]
longitude = ds["longitude"][:]
z = ds["Z"][:]

skip = 1
var = field[begin:skip:end, begin:skip:end, :, 1]
# lat-lon
lat = latitude[begin:skip:end, 1]
lon = longitude[begin:skip:end, 1]

# Set missing values to NaN
missing_value = 0.0
bools = (var .== missing_value) .| (var .=== NaN) .| (var .=== missing)
var[bools] .= NaN

function available_fields(; local_directory="/mnt/podaac_drive/Version4/Release4/interp_monthly/")
    afields = readdir(local_directory)
    println("The availabel fields are ", afields)
    return afields
end

function grab_ecco_dataset(fieldname::String, year_int::Int, month_int::Int;
    local_directory="/mnt/podaac_drive/Version4/Release4/interp_monthly/")

    fields = readdir(local_directory)
    if !(fieldname in fields)
        println(fieldname, " is not an available field")
        println("please use ", fields)
        return fields
    end
    # Find available years
    years = readdir(local_directory * fieldname)
    year = years[year_int]
    # Find available months
    println("Looking at year ", year)
    months = readdir(joinpath(local_directory * fieldname, year))
    month = months[month_int]
    tic = time()
    println("Looking at month file ", month)
    full_path = joinpath(local_directory, fieldname, year, month)
    # Load Data Set
    ds = Dataset(full_path, "r")
    toc = time()
    println("Time to load dataset = ", toc - tic, " seconds")

    return ds
end

function grab_ecco_field(fieldname::String,  year_int::Int, month_int::Int;
    skip=1, local_directory="/mnt/podaac_drive/Version4/Release4/interp_monthly/")

    ds = grab_ecco_dataset(fieldname, year_int, month_int, local_directory=local_directory)

    field = ds[fieldname]

    # Plot variables
    latitude = ds["latitude"][:]
    longitude = ds["longitude"][:]
    z = ds["Z"][:]

    skip = 1
    var = field[begin:skip:end, begin:skip:end, :, 1]
    # lat-lon
    lat = latitude[begin:skip:end, 1]
    lon = longitude[begin:skip:end, 1]

    # Set missing values to NaN
    missing_value = 0.0
    bools = (var .== missing_value) .| (var .=== NaN) .| (var .=== missing)
    var[bools] .= NaN

    return var, lat, lon, z
end