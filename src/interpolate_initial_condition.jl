function interp_initial_condition(old_vector, Nx_old, Ny_old, Nx_new, Ny_new, Nz, oldlat, lat)

    old_grid = RectilinearGrid(size = (Nx_old, Ny_old, 1), y = (-oldlat, oldlat), x = (-180, 180), z = (0, 1))
    new_grid = RectilinearGrid(size = (Nx_new, Ny_new, 1), y = (-lat, lat), x = (-180, 180), z = (0, 1))

    new_vector = zeros(Nx_new, Ny_new, Nz)
    old_vec = zeros(Nx_old, Ny_old)
    for k = 1:Nz
        old_vec = old_vector[:, :, k]
        new_vec = interpolate_one_level(old_vec, old_grid, new_grid, Center)
        new_vector[:, :, k] = new_vec
    end

    return new_vector
end

function check_zeros(bathymetry, z_faces, old_array; ensure_positivity = false, max_passes = 10)
    Nx, Ny, Nz = size(old_array)
    grid = RectilinearGrid(size = (Nx, Ny, Nz), x = (0, 1), y = (0, 1), topology = (Periodic, Periodic, Bounded), z = z_faces, halo = (1, 1, 1))

    zc = grid.zᵃᵃᶜ[1:Nz]

    if ensure_positivity 
        old_array[old_array .< 0] .= 0.0
    end

    field = Field{Center, Center, Center}(grid)
    set!(field, old_array)
    fill_halo_regions!(field, architecture(grid))
    
    condition = false
    pass = 1
    while condition == false && pass < max_passes + 1
        @info "Pass number $pass"
        condition = true    
        for k in 1:Nz
            @info "we are at k = $k"
            for i in 1:Nx, j in 1:Ny
                if zc[k] > bathymetry[i, j]
                    if old_array[i, j, k] == 0
                        condition = false
                        @info "There is a zero! at $i, $j with z = $(zc[k]) and bat = $(bathymetry[i, j])"
                        neigbours = [field[i + 1, j, k], field[i - 1, j, k], field[i, j + 1, k], field[i, j - 1, k], field[i, j, k + 1]]
                        non_null  = Int.(neigbours .!= 0)
                        @info "Old array: $(old_array[i, j, k]), sum neigbours $(sum(neigbours))"
                        if sum(non_null) != 0 
                            old_array[i, j, k] = dot(non_null, neigbours) / sum(non_null)
                        end
                    end
                end
            end
        end
        pass += 1
    end

    return old_array
end

function substitute_min_value(old_array, min_val) 
    old_array[old_array.< min_val] .= min_val 
    return old_array
end


function read_and_interpolate_var(var, folder, Nxₒ, Nyₒ, Nxₙ, Nyₙ, Nzₙ, latitude)

    println("$var")
            
    path_to_file = "/mnt/podaac_drive/Version4/Release4/interp_monthly/$folder/1992/$(folder)_1992_01.nc"

    interp_full = ncread(path_to_file, folder)

    interp_full = reverse(interp_full[:, :, :, 1], dims = 3)
    
    passes        = 20
    filter_passes = 10

    interp_old = deepcopy(interp_full)

    Nx, Ny, Nz = size(interp_full)

    for pass = 1:passes
        @info "pass $pass in propagate step"
        interp_full = propagate_step(interp_full, Nz)
    end

    for pass in 1:filter_passes
        @info "pass $pass in horizontal filter"
        interp_full = horizonthal_filter(interp_full, 2, Nx-1, 2, Ny-1, Nz)
    end

    # # We have to homogenize boundary values which seem quite off (periodic boundaries)
    # interp_full = filter_boundary_values(interp_full, Nz, passes) 

    # Remove boundary values
    edges = [collect(1:3)..., [Nx - i for i = 0:2]...]
    subs  = [[4 for i in 1:3]..., [Nx-3 for i in 1:3]...]
    
    for (edge, sub) in zip(edges, subs)
        interp_full[edge, :, :] = interp_full[sub, :, :] 
    end

    # We have to homogenize boundary values which seem quite off (periodic boundaries)
    interp_full = filter_boundary_values(interp_full, Nz, filter_passes) 

    interp_new  = interp_initial_condition(interp_full, Nxₒ, Nyₒ, Nxₙ, Nyₙ, Nzₙ, old_latitude, latitude)


    # Remove boundary values
    edges = [collect(1:3)..., [Nxₙ - i for i = 0:2]...]
    subs  = [[4 for i in 1:3]..., [Nxₙ-3 for i in 1:3]...]
    
    for (edge, sub) in zip(edges, subs)
        interp_new[edge, :, :] = interp_new[sub, :, :] 
    end
    
    # We have to homogenize boundary values which seem quite off (periodic boundaries)
    interp_new = filter_boundary_values(interp_new, Nz, 5) 

    return (interp_new, interp_old)
end

function read_and_interpolate_var_from_array(old_array, Nxₙ, Nyₙ, Nzₙ, old_latitude, latitude)

    passes        = 20
    filter_passes = 10

    interp_old = deepcopy(old_array)

    Nx, Ny, Nz = size(interp_old)

    for pass = 1:passes
        @info "pass $pass in propagate step"
        interp_old = propagate_step(interp_old, Nz)
    end

    for pass in 1:filter_passes
        @info "pass $pass in horizontal filter"
        interp_old = horizonthal_filter(interp_old, 2, Nx-1, 2, Ny-1, Nz)
    end

    # Remove boundary values
    edges = [collect(1:3)..., [Nx - i for i = 0:2]...]
    subs  = [[4 for i in 1:3]..., [Nx-3 for i in 1:3]...]
    
    for (edge, sub) in zip(edges, subs)
        interp_old[edge, :, :] = interp_old[sub, :, :] 
    end

    # We have to homogenize boundary values which seem quite off (periodic boundaries)
    interp_old    = filter_boundary_values(interp_old, Nz, filter_passes) 
    interp_array  = interp_initial_condition(interp_old, Nx, Ny, Nxₙ, Nyₙ, Nzₙ, old_latitude, latitude)

    # Remove boundary values
    edges = [collect(1:3)..., [Nxₙ - i for i = 0:2]...]
    subs  = [[4 for i in 1:3]..., [Nxₙ-3 for i in 1:3]...]
    
    for (edge, sub) in zip(edges, subs)
        interp_array[edge, :, :] = interp_array[sub, :, :] 
    end
    
    # We have to homogenize boundary values which seem quite off (periodic boundaries)
    interp_array = filter_boundary_values(interp_array, Nz, 5) 

    return interp_array
end

# vars    = ("temp", "salt")
# folders = ("THETA", "SALT")

# for (var, folder) in zip(vars, folders)

#     latitude = 75
#     degree = 1/12

#     Nxₒ  = 1440
#     Nyₒ  = 600

#     Nxₙ  = Int(360 / degree)
#     Nyₙ  = Int(2latitude / degree)
    
#     Nzₒ  =  Nzₙ  = 50

#     interpolated, smoothened, original = read_and_interpolate_var(var, folder, Nxₒ, Nyₒ, Nxₙ, Nyₙ, Nzₙ, latitude)

#     bathymetry   = jldopen("../bathymetry/bathymetry-$(Nxₙ)x$(Nyₙ)-latitude-$(latitude).jld2")["bathymetry"]
#     z_faces      = jldopen("../bathymetry/z_faces-$(Nzₙ)-levels.jld2")["z_faces"]
    
#     # @info "Checking that there are no zeros at z > bathymetry"
#     interpolated = check_zeros(bathymetry, z_faces, interpolated)

#     output_i = var * "_interpolated_$(Nxₙ)x$(Nyₙ)x$(Nzₙ)-latitude-$(latitude)"
#     output_s = var * "_smoothened_$(Nxₒ)x$(Nyₒ)x$(Nzₙ)-latitude-90"
#     output_o = var * "_original_$(Nxₒ)x$(Nyₒ)x$(Nzₙ)-latitude-90"
#     jldsave(output_i * ".jld2", field = interpolated)
#     jldsave(output_s * ".jld2", field = smoothened)
#     jldsave(output_o * ".jld2", field = original)
#     field_to_bin = bswap.(interpolated)
    
#     write(output_i * ".bin", field_to_bin)
# end
