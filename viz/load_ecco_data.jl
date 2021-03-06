using HDF5, Statistics

function load_ecco_data(filepath, month_index)
    fid = h5open(filepath, "r")

    ϕ = read(fid["grid"]["latitude"])
    λ = read(fid["grid"]["longitude"])
    z = read(fid["grid"]["raditude"])

    ecco_salt = read(fid["fields"]["salinity"])
    ecco_temperature = read(fid["fields"]["temperature"])

    ecco_S = ecco_salt[:, :, :, month_index]
    ecco_θ = ecco_temperature[:, :, :, month_index]
    ecco_grid = (; λ, ϕ, z)
    close(fid)

    return ecco_grid, ecco_S, ecco_θ
end

# seasonal cycle statistics
# assumed shape [ix, iy, iz, iω]
function return_statistics(field, filepath)
    load_ecco_data(filepath)
    ecco_size = size(field)
    r_field = reshape(field, (prod(ecco_size[1:3]), ecco_size[4]))
    upper_quantile = zeros(ecco_size[1:3])
    lower_quantile = zeros(ecco_size[1:3])
    median_quantile = zeros(ecco_size[1:3])
    for i in 1:prod(ecco_size[1:3])
        if isnan(r_field[i, 1])
            mean_salt[i] = NaN
            upper_quantile[i] = NaN
            lower_quantile[i] = NaN
        else
            upper_quantile[i] = quantile(r_field[i, :], 0.75)
            lower_quantile[i] = quantile(r_field[i, :], 0.25)
            median_quantile[i] = quantile(r_field[i, :], 0.5) # perhaps have nice_salt?
        end
    end
    return upper_quantile, median_quantile, lower_quantile 
end

