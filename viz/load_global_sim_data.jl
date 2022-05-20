using Printf, Oceananigans, JLD2, SeawaterPolynomials

latitude = (-75, 75)

function load_global_sim_data(; global_filepath = "/storage1/uq-global/GlobalShenanigans.jl/", 
                                load = :grid, size = (1, 1, 1),
                                ic_filepath = "longer_null_hypothesis_teos10.jld2")
    
    Nx, Ny, Nz = size
    if load == :simulation
        println("loading from simulation")
        uvel = Array(interior(simulation.model.velocities.u))
        vvel = Array(interior(simulation.model.velocities.v))
        Ttrac = Array(interior(simulation.model.tracers.T))
        Strac = Array(interior(simulation.model.tracers.S))
    else
        # 1 degree resolution
        arch = CPU()
        # A spherical domain
        z_faces = jldopen("data/zgrid.jld2")["z"][5:end-4]
        grid = LatitudeLongitudeGrid(arch,
            size=(Nx, Ny, Nz),
            longitude=(-180, 180),
            latitude=latitude,
            halo=(4, 4, 4),
            z=z_faces,
            precompute_metrics=true)
        if load == :file
            println("loading from ", ic_filepath)
            jlfile = jldopen(ic_filepath)
            uvel = jlfile["velocities"]["u"]
            vvel = jlfile["velocities"]["v"]
            Ttrac = jlfile["tracers"]["T"]
            Strac = jlfile["tracers"]["S"]
            eta_fs = jlfile["free_surface"]["eta"]
            close(jlfile)
        end
    end

    z = Array(grid.zᵃᵃᶜ[1:end-4])

    eos_teos10 = SeawaterPolynomials.TEOS10EquationOfState()
    rz = reshape(z, (1, 1, length(z)))

    if !(load == :grid)
        # 2km reference level for potential density
        teos10_ρ = SeawaterPolynomials.ρ′.(Ttrac, Strac, Ref(-2000.0), Ref(eos_teos10))
        ρ_prime = @. (teos10_ρ - eos_teos10.reference_density) / eos_teos10.reference_density

        ρ_min, ρ_max = extrema(ρ_prime[(!).(isnan.(ρ_prime))])
        levels = range(ρ_min, ρ_max, length=41)

        continentsu = (uvel .== 0.0)
        continentsv = (vvel .== 0.0)
        continentst = (Ttrac .== 0.0)
        uvel[continentsu] .= NaN
        vvel[continentsv] .= NaN
        Ttrac[continentst] .= NaN
        Strac[continentst] .= NaN
        ρ_prime[continentst] .= NaN
    end

    z = Array(grid.zᵃᵃᶜ[1:end-4])
    λc = Array(grid.λᶜᵃᵃ[1:end-4])
    φc = Array(grid.φᵃᶜᵃ[1:end-4])
    φf = Array(grid.φᵃᶠᵃ[1:end-4])

    oceananigans_grid = (; λ = λc, ϕ = φc, z = z)

    return oceananigans_grid
end