using Oceananigans
using Oceananigans.Operators: ζ₃ᶠᶠᶜ
using JLD2

function compute_vorticity!(folder, size)

    Nx, Ny, Nz = size
    latitude = (-75, 75)

    output_prefix = folder * "near_global_lat_lon_$(Nx)_$(Ny)_$(Nz)_fine"

    z_faces = jldopen("data/z_faces-50-levels.jld2")["z_faces"]
    z_faces = z_faces[3:end]

    # A spherical domain
    @show grid = LatitudeLongitudeGrid(CPU(),
                                    size = (Nx, Ny, Nz),
                                    longitude = (-180, 180),
                                    latitude = latitude,
                                    halo = (3, 3, 3),
                                    z = z_faces,
                                    precompute_metrics = true)

    surface_file = jldopen(output_prefix * "_surface.jld2")
    iterations   = parse.(Int, keys(surface_file["timeseries/t"]))

    ζ₃ = zeros(Nx, Ny+1, length(iterations))
    u  = Field{Face, Center, Nothing}(grid)
    v  = Field{Center, Face, Nothing}(grid)

    for (time, iter) in enumerate(iterations)
    uvec = surface_file["timeseries/u/$iter"][:, :, end]
        vvec = surface_file["timeseries/v/$iter"][:, :, end]
        Oceananigans.set!(u, uvec)
        Oceananigans.set!(v, vvec)
        
        for i in 1:Nx, j in 1:Ny+1
            ζ₃[i, j, time] = ζ₃ᶠᶠᶜ(i, j, 1, grid, u, v)
        end
    end

    jldsave(folder * "vorticity.jld2", vort = ζ₃)
end

function compute_vorticity!(folder, size_grid, iteration)

    Nx, Ny, Nz = size_grid
    latitude = (-75, 75)

    output_prefix = folder * "near_global_lat_lon_$(Nx)_$(Ny)_$(Nz)_fine_checkpoint_iteration$iteration"

    z_faces = jldopen("data/z_faces-50-levels.jld2")["z_faces"]
    z_faces = z_faces[3:end]

    # A spherical domain
    @show grid = LatitudeLongitudeGrid(CPU(),
                                    size = (Nx, Ny, Nz),
                                    longitude = (-180, 180),
                                    latitude = latitude,
                                    halo = (3, 3, 3),
                                    z = z_faces,
                                    precompute_metrics = true)

    velocities = jldopen(output_prefix * ".jld2")

    ζ₃ = zeros(Nx, Ny+1, Nz)
    u  = Field{Face, Center, Nothing}(grid)
    v  = Field{Center, Face, Nothing}(grid)

    uvec = velocities["u/data"][5:end-4, 5:end-4, end-4]
    vvec = velocities["v/data"][5:end-4, 5:end-4, end-4]

    Oceananigans.set!(u, uvec)
    Oceananigans.set!(v, vvec)

    ke = Field(u^2 + v^2)
    compute!(ke)

    δ  = Field(∂x(u) + ∂y(v))
    compute!(δ)
     
    div = interior(δ)
    for i = 1:Nx, j = 1:Ny+1 
        ζ₃[i, j, 1] = ζ₃ᶠᶠᶜ(i, j, Nz, grid, u, v)
    end

    jldsave(folder * "vorticity$iteration.jld2", vort = ζ₃, div = div, ke = ke)
end
