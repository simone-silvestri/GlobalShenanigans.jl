# grab interpolate_fluxes
function read_and_interpolate_flux(var, folder, yr, dims)
    lat = 75
    degree = 1/4

    Nxₒ  = 720
    Nyₒ  = 360
    Nxₙ  = Int(360 / degree)
    Nyₙ  = Int(2lat/ degree)

    old   = 90
    gridₒ = RectilinearGrid(size = (Nxₒ, Nyₒ, 1), y = (-old, old) , x = (-180, 180), z = (0, 1))
    gridₙ = RectilinearGrid(size = (Nxₙ, Nyₙ, 1), y = (-lat, lat), x = (-180, 180), z = (0, 1))

    xₒ = gridₒ.xᶜᵃᵃ[1:Nxₒ]
    yₒ = gridₒ.yᵃᶜᵃ[1:Nyₒ]
    xₙ = gridₙ.xᶜᵃᵃ[1:Nxₙ]
    yₙ = gridₙ.yᵃᶜᵃ[1:Nyₙ]
    
    flux = zeros(Nxₙ, Nyₙ, 12)
    months  = ("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

    for (idx, month) in enumerate(months)
        println("$var, $yr, $month")
        
        path_to_file = "/mnt/podaac_drive/Version4/Release4/interp_monthly/$folder/$yr/$(folder)_$(yr)_$(month).nc"

        interp = ncread(path_to_file, folder)
        for pass = 1:10
            @info "pass $pass of propagate step"
            interp = propagate_step(interp, 1)
        end

        if dims ==3
            for pass = 1:100
                @info "pass $pass of propagate step"
                interp = propagate_step(interp, 1)
            end
            interp = interp[:, :, :, 1]
            interp = filter_boundary_values(interp, 1, 10) 
            interp = interp[:, :, 1]
        end

        interp_full = deepcopy(interp)
        interp = interpolate_one_level(interp_full, gridₒ, gridₙ, Center) 

        flux[:, :, idx] .= interp[:, :]
    end

    @show size(flux)
    return flux
end



for years in collect(2001:2007)
    vars    = ("tau_x", "tau_y")
    folders = ("oceTAUE", "oceTAUN")

    for (var, folder) in zip(vars, folders)

        interp = read_and_interpolate_flux(var, folder, years, 2)

        jldsave(var * "-1440x600-latitude-75-$years.jld2", field = interp)
        field_to_bin = bswap.(interp)
        
        write(var * "-1440x600-latitude-75-$years.bin", field_to_bin)
    end

    vars = ("temp", "salt")
    folders = ("THETA", "SALT")


    for (var, folder) in zip(vars, folders)
        interp = read_and_interpolate_flux(var, folder, years, 3)

        Nx = size(interp, 1)

        # Remove boundary values
        edges = [collect(1:3)..., [Nx - i for i = 0:2]...]
        subs  = [[4 for i in 1:3]..., [Nx-3 for i in 1:3]...]
        
        for (edge, sub) in zip(edges, subs)
            interp[edge, :, :] = interp[sub, :, :] 
        end

        jldsave(var * "-1440x600-latitude-75-$years.jld2", field = interp)
        field_to_bin = bswap.(interp)
        
        write(var * "-1440x600-latitude-75-$years.bin", field_to_bin)
    end
end
