using HDF5
fname = "wrangle_ecco_2.h5"
fid = h5open(fname, "r")

ϕ = read(fid["grid"]["latitude"])
λ = read(fid["grid"]["longitude"])
z = read(fid["grid"]["raditude"])

ecco_salt = read(fid["fields"]["salinity"])
ecco_temperature = read(fid["fields"]["temperature"])

month_index = 1

ecco_S = view(ecco_salt, :, :, :, month_index)
ecco_θ = view(ecco_temperature, :, :, :, month_index)
ecco_grid = (; λ, ϕ, z)
close(fid)