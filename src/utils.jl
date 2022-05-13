@inline function z_face_from_center(z_center)

    z_faces = zeros(length(z_center) + 1)
    for k = 2:length(z_center) + 1
        z_faces[k] = z_faces[k-1] + 2*(z_center[k-1] - z_faces[k-1])
    end

    return reverse(z_faces)
end

@inline function interpolate_one_level(old_array, old_grid, new_grid, loc)
    old_field = Field{Center, Center, Center}(old_grid)
    set!(old_field, old_array)

    Nx_new, Ny_new = size(new_grid)[[1, 2]]
    new_array = zeros(Nx_new, Ny_new)

    for i in 1:Nx_new, j in 1:Ny_new
        if loc == Face
            new_array[i, j] = interpolate(old_field, new_grid.xᶜᵃᵃ[i], new_grid.yᵃᶠᵃ[j], old_grid.zᵃᵃᶜ[1])
        else
            new_array[i, j] = interpolate(old_field, new_grid.xᶜᵃᵃ[i], new_grid.yᵃᶜᵃ[j], old_grid.zᵃᵃᶜ[1])
        end
    end
    return new_array
end

@inline function propagate_step(vec, Nz_mine)
    Nx, Ny, Nz = size(vec)
    vec2 = deepcopy(vec)
    for k in 1:Nz_mine
        for i in 2:Nx-1, j in 2:Ny-1
            neigbours = [vec[i + 1, j, k], vec[i - 1, j, k], vec[i, j + 1, k], vec[i, j - 1, k]]
            non_null  = Int.(neigbours .!= 0)
            if (vec[i, j, k] == 0) && sum(non_null) > 0
                vec2[i, j, k] = dot(non_null, neigbours) / sum(non_null) 
            end
        end
    end
    return vec2
end

@inline function horizonthal_filter(vec, iᵢ, iₑ, jᵢ, jₑ, Nz_mine)
    vec2 = deepcopy(vec)
    for k in 1:Nz_mine
        for i in iᵢ:iₑ, j in jᵢ:jₑ
            neigbours = [vec[i, j, k], vec[i + 1, j, k], vec[i - 1, j, k], vec[i, j + 1, k], vec[i, j - 1, k]]
            non_null  = Int.(neigbours .!= 0)
            if sum(non_null) > 0
                vec2[i, j, k] = sum(neigbours) / 5
            end
        end
    end
    return vec2
end

@inline function filter_boundary_values(old_vector, Nz_mine, passes)

    Nx, Ny, Nz = size(old_vector)
    new_vector = deepcopy(old_vector)
    grid = RectilinearGrid(size = (Nx, Ny), 
                              y = (-1, 1),
                              x = (-1, 1),
                       topology = (Periodic, Bounded, Flat), halo = (20, 20))

    field = Field{Center, Center, Nothing}(grid)
    for k = 1:Nz_mine
        set!(field, old_vector[:, :, k])
        fill_halo_regions!(field, architecture(grid))
        @info "level $k in filter boundary"
        for pass = 1:passes
            field = horizonthal_filter(field,   -15,    30, 2, Ny-1, 1)
            field = horizonthal_filter(field, Nx-30, Nx+15, 2, Ny-1, 1)

            new_vector[:, :, k] = interior(field)[:, :, 1]
        end
    end
    return new_vector
end

@inline function copy_boundary_values(old_vector, Nz_mine)

    Nx, Ny, Nz = size(old_vector)
    new_vector = deepcopy(old_vector)
    grid = RectilinearGrid(size = (Nx, Ny), 
                              y = (-1, 1),
                              x = (-1, 1),
                       topology = (Periodic, Bounded, Flat), halo = (20, 20))

    field = Field{Center, Center, Nothing}(grid)
    for k = 1:Nz_mine
        set!(field, old_vector[:, :, k])
        fill_halo_regions!(field, architecture(grid))
        @info "level $k in filter boundary"
        for pass = 1:100
            field = horizonthal_filter(field,   -15,    30, 2, Ny-1, 1)
            field = horizonthal_filter(field, Nx-30, Nx+15, 2, Ny-1, 1)

            new_vector[:, :, k] = interior(field)[:, :, 1]
        end
    end
    return new_vector
end