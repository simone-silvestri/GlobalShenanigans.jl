using GlobalShenanigans, JLD2
using GLMakie

function visualize_bathymetry(bathy)
    bat = deepcopy(bathy)
    bat[ bat .== 100 ] .= NaN
    heatmap(bat, nan_colors = :black)
end

degree   = 1.0 
latitude = 75.0

Nx = Int(360 / degree)
Ny = Int(2latitude / degree)
Nz = 48

passes        = 10
length_lake_x = 20
length_lake_y = 20

## Start of the interpolation
arch = CPU()
bathymetry_original = GlobalShenanigans.interpolate_bathymetry_from_file("data/bathymetry-ice-21600x10800.jld2", passes, degree, latitude)

## Start the modifications
bathy = Vector(undef, 5)

# mask individual elements
# bathy[1] = GlobalShenanigans.mask_isolated_elements(bathymetry_original, arch)
# bathy[2] = GlobalShenanigans.mask_coupled_elements(bathy[1], arch)
# bathy[3] = GlobalShenanigans.fill_up_all_lakes(bathy[2], arch, length_lake_x, length_lake_y)

# remove caspian sea
bathy[1][225:240, 110:125] .= ABOVE_SEA_LEVEL

# remove baltic sea
bathy[1][198:210, 129:141] .= ABOVE_SEA_LEVEL
bathy[1][194:198, 129:132] .= ABOVE_SEA_LEVEL
bathy[1][214:221, 139:142] .= ABOVE_SEA_LEVEL

# remove black sea
bathy[1][210:225, 117:122] .= ABOVE_SEA_LEVEL

# opening up mediterranean
bathy[1][175, 111:112] .= 0.5 * (bathy[1][174, 111:112] .+ bathy[1][176, 111:112])

# opening up red sea
bathy[1][222:223, 89] .= bathy[1][222:223, 90]
bathy[1][224, 89]      = bathy[1][224, 88]

GlobalShenanigans.write_bathymetry_to_file("bathymetry", bathy[1], latitude)

using DataDeps

path = "https://github.com/CliMA/OceananigansArtifacts.jl/raw/ss/new_hydrostatic_data_after_cleared_bugs/quarter_degree_near_global_input_data/"

datanames = ["temp-1440x600-latitude-75",
             "salt-1440x600-latitude-75",
             "tau_x-1440x600-latitude-75",
             "tau_y-1440x600-latitude-75",
             "initial_conditions"]

dh = DataDep("quarter_degree_near_global_lat_lon",
    "Forcing data for global latitude longitude simulation",
    [path * data * ".jld2" for data in datanames]
)

DataDeps.register(dh)

datadep"quarter_degree_near_global_lat_lon"

files = [:file_temp, :file_salt, :file_tau_x, :file_tau_y, :file_init]
for (data, file) in zip(datanames, files)
    datadep_path = @datadep_str "quarter_degree_near_global_lat_lon/" * data * ".jld2"
    @eval $file = jldopen($datadep_path)
end

# T = file_init["T"]
# S = file_init["S"]

# T_interp = GlobalShenanigans.read_and_interpolate_var_from_array(T, Nx, Ny, Nz, latitude, latitude)
# S_interp = GlobalShenanigans.read_and_interpolate_var_from_array(S, Nx, Ny, Nz, latitude, latitude)

# jldsave("initial_conditions-1degree.jld2", T = T_interp, S = S_interp)

τˣ_old = file_tau_x["field"]
τʸ_old = file_tau_y["field"]
Tˢ_old = file_temp["field"]
Sˢ_old = file_salt["field"]
τˣ = zeros(Nx, Ny, 12) 
τʸ = zeros(Nx, Ny, 12)
Tˢ = zeros(Nx, Ny, 12) 
Sˢ = zeros(Nx, Ny, 12) 

τˣ = GlobalShenanigans.read_and_interpolate_var_from_array(τˣ_old, Nx, Ny, 12, latitude, latitude)
τʸ = GlobalShenanigans.read_and_interpolate_var_from_array(τʸ_old, Nx, Ny, 12, latitude, latitude)
Tˢ = GlobalShenanigans.read_and_interpolate_var_from_array(Tˢ_old, Nx, Ny, 12, latitude, latitude)
Sˢ = GlobalShenanigans.read_and_interpolate_var_from_array(Sˢ_old, Nx, Ny, 12, latitude, latitude)

jldsave("boudary_conditions-1degree.jld2", τˣ = τˣ, τʸ = τʸ, Tˢ = Tˢ, Sˢ = Sˢ)
