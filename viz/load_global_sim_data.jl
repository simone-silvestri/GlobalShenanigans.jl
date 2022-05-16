using GLMakie, Printf, Oceananigans, JLD2, SeawaterPolynomials

load_from_simulation = false
if load_from_simulation
    println("loading from simulation")
    uvel = Array(interior(simulation.model.velocities.u))
    vvel = Array(interior(simulation.model.velocities.v))
    Ttrac = Array(interior(simulation.model.tracers.T))
    Strac = Array(interior(simulation.model.tracers.S))
else
    ic_filepath = "longer_null_hypothesis_teos10.jld2"
    println("loading from ", ic_filepath)
    jlfile = jldopen(ic_filepath)
    uvel = jlfile["velocities"]["u"]
    vvel = jlfile["velocities"]["v"]
    Ttrac = jlfile["tracers"]["T"]
    Strac = jlfile["tracers"]["S"]
    eta_fs = jlfile["free_surface"]["eta"]
    close(jlfile)

    z_faces = jldopen("zgrid.jld2")["z"][5:end-4]

    latitude = (-75, 75)

    # 1 degree resolution
    Nx = 360
    Ny = 150
    Nz = 48
    arch = CPU()
    # A spherical domain
    z_faces = jldopen("zgrid.jld2")["z"][5:end-4]
    grid = LatitudeLongitudeGrid(arch,
        size=(Nx, Ny, Nz),
        longitude=(-180, 180),
        latitude=latitude,
        halo=(4, 4, 4),
        z=z_faces,
        precompute_metrics=true)

end

z = Array(grid.zᵃᵃᶜ[1:end-4])
eos_teos10 = SeawaterPolynomials.TEOS10EquationOfState()
rz = reshape(z, (1,1,length(z)))
# 2km reference level for potential density
teos10_ρ = SeawaterPolynomials.ρ′.(Ttrac, Strac, Ref(-2000.0), Ref(eos_teos10))
ρ_prime = @. (teos10_ρ - eos_teos10.reference_density) / eos_teos10.reference_density

ρ_min, ρ_max = extrema(ρ_prime[(!).(isnan.(ρ_prime))])
levels = range(ρ_min, ρ_max, length=41)

z = Array(grid.zᵃᵃᶜ[1:end-4])
λc = Array(grid.λᶜᵃᵃ[1:end-4])
φc = Array(grid.φᵃᶜᵃ[1:end-4])
φf = Array(grid.φᵃᶠᵃ[1:end-4])

continentsu = (uvel .== 0.0)
continentsv = (vvel .== 0.0)
continentst = (Ttrac .== 0.0)
uvel[continentsu] .= NaN
vvel[continentsv] .= NaN
Ttrac[continentst] .= NaN
Strac[continentst] .= NaN
ρ_prime[continentst] .= NaN
