#import Pkg
#Pkg.activate("/home/ssilvest/Oceananigans.jl/")
using Statistics
using JLD2
using Printf
# using Plots
using Oceananigans
using Oceananigans.Units

using Oceananigans.Fields: interpolate, Field
using Oceananigans.Architectures: arch_array
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis
using Oceananigans.BoundaryConditions
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom, solid_node, solid_interface
using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity
using CUDA: @allowscalar, device!
using Oceananigans.Operators
using Oceananigans.Operators: Δzᵃᵃᶜ
using Oceananigans: prognostic_fields
using Oceananigans.Advection: EnergyConservingScheme, VelocityStencil, VorticityStencil

device!(0)

@inline function visualize(field, lev, dims) 
    (dims == 1) && (idx = (lev, :, :))
    (dims == 2) && (idx = (:, lev, :))
    (dims == 3) && (idx = (:, :, lev))

    r = deepcopy(Array(interior(field)))[idx...]
    r[ r.==0 ] .= NaN
    return r
end

#####
##### Grid
#####

arch = GPU()
reference_density = 1029

latitude = (-75, 75)

# 1 degree resolution
Nx = 360
Ny = 150
Nz = 48

const Nyears  = 5
const Nmonths = 12 
const thirty_days = 30days

output_prefix = "near_global_lat_lon_$(Nx)_$(Ny)_$(Nz)_fine"
pickup_file   = false 

#####
##### Load forcing files and inital conditions from ECCO version 4
##### https://ecco.jpl.nasa.gov/drive/files
##### Bathymetry is interpolated from ETOPO1 https://www.ngdc.noaa.gov/mgg/global/
#####
#=
using DataDeps

path = "https://github.com/CliMA/OceananigansArtifacts.jl/raw/ss/new_hydrostatic_data_after_cleared_bugs/quarter_degree_near_global_input_data/"

datanames = "z_faces-50-levels"

dh = DataDep("quarter_degree_near_global_lat_lon",
    "Forcing data for global latitude longitude simulation",
    [path * data * ".jld2" for data in datanames]
)

DataDeps.register(dh)

datadep"quarter_degree_near_global_lat_lon"

datadep_path = @datadep_str "quarter_degree_near_global_lat_lon/z_faces-50-levels.jld2"
=#

bathymetry = jldopen("bathymetry-360x150-latitude-75.0.jld2")["bathymetry"]

τˣ = zeros(Nx, Ny, Nmonths)
τʸ = zeros(Nx, Ny, Nmonths)
T★ = zeros(Nx, Ny, Nmonths)
S★ = zeros(Nx, Ny, Nmonths)

# Files contain 1 year (1992) of 12 monthly averages
τˣ = jldopen("boundary_conditions-1degree.jld2")["τˣ"] ./ reference_density
τʸ = jldopen("boundary_conditions-1degree.jld2")["τˣ"] ./ reference_density
T★ = jldopen("boundary_conditions-1degree.jld2")["Tˢ"] 
S★ = jldopen("boundary_conditions-1degree.jld2")["Sˢ"] 

# Remember the convention!! On the surface a negative flux increases a positive decreases
bathymetry = arch_array(arch, bathymetry)
τˣ = arch_array(arch, - τˣ)
τʸ = arch_array(arch, - τʸ)

target_sea_surface_temperature = T★ = arch_array(arch, T★)
target_sea_surface_salinity    = S★ = arch_array(arch, S★)

# Stretched faces taken from ECCO Version 4 (50 levels in the vertical)
z_faces = jldopen("zgrid.jld2")["z"][5:end-4]

# A spherical domain
@show underlying_grid = LatitudeLongitudeGrid(arch,
                                              size = (Nx, Ny, Nz),
                                              longitude = (-180, 180),
                                              latitude = latitude,
                                              halo = (4, 4, 4),
                                              z = z_faces,
                                              precompute_metrics = true)

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry))

#####
##### Physics and model setup
#####

νh = 1e+1 * 1e5
νz = 5e-3 * 1e2
κh = 1e+1 * 1e5
κz = 1e-4 * 1e2

using Oceananigans.Operators: Δx, Δy
using Oceananigans.TurbulenceClosures

@inline νhb(i, j, k, grid, lx, ly, lz) = (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2 ))^2 / 5days

horizontal_diffusivity = HorizontalScalarDiffusivity(ν=νh, κ=κh)
vertical_diffusivity   = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=νz, κ=κz)
convective_adjustment  = RiBasedVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), κ₀ = 1.0)
biharmonic_viscosity   = HorizontalScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true)
                                                    
#####
##### Boundary conditions / time-dependent fluxes 
#####

@inline current_time_index(time, tot_months)     = mod(unsafe_trunc(Int32, time / thirty_days), tot_months) + 1
@inline next_time_index(time, tot_months)        = mod(unsafe_trunc(Int32, time / thirty_days) + 1, tot_months) + 1
@inline cyclic_interpolate(u₁::Number, u₂, time) = u₁ + mod(time / thirty_days, 1) * (u₂ - u₁)

Δz_top    = @allowscalar Δzᵃᵃᶜ(1, 1, grid.Nz, grid.grid)
Δz_bottom = @allowscalar Δzᵃᵃᶜ(1, 1, 1, grid.grid)

@inline function surface_wind_stress(i, j, grid, clock, fields, τ)
    time = clock.time
    n₁ = current_time_index(time, Nmonths)
    n₂ = next_time_index(time, Nmonths)

    @inbounds begin
        τ₁ = τ[i, j, n₁]
        τ₂ = τ[i, j, n₂]
    end

    return cyclic_interpolate(τ₁, τ₂, time)
end

u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τˣ)
v_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τʸ)

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1]
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1]

# Linear bottom drag:
μ = 0.001 # ms⁻¹

@inline is_immersed_drag_u(i, j, k, grid) = Int(solid_interface(Face(), Center(), Center(), i, j, k-1, grid) & !solid_node(Face(), Center(), Center(), i, j, k, grid))
@inline is_immersed_drag_v(i, j, k, grid) = Int(solid_interface(Center(), Face(), Center(), i, j, k-1, grid) & !solid_node(Center(), Face(), Center(), i, j, k, grid))                                

# Keep a constant linear drag parameter independent on vertical level
@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * is_immersed_drag_u(i, j, k, grid) * fields.u[i, j, k] / Δzᵃᵃᶜ(i, j, k, grid)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * is_immersed_drag_v(i, j, k, grid) * fields.v[i, j, k] / Δzᵃᵃᶜ(i, j, k, grid)

Fu = Forcing(u_immersed_bottom_drag, discrete_form = true, parameters = μ)
Fv = Forcing(v_immersed_bottom_drag, discrete_form = true, parameters = μ)

u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ)
v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ)

@inline function surface_temperature_relaxation(i, j, grid, clock, fields, p)
    time = clock.time

    n₁ = current_time_index(time, Nmonths)
    n₂ = next_time_index(time, Nmonths)

    @inbounds begin
        T★₁ = p.T★[i, j, n₁]
        T★₂ = p.T★[i, j, n₂]
        T_surface = fields.T[i, j, grid.Nz]
    end

    T★ = cyclic_interpolate(T★₁, T★₂, time)
                                
    return p.λ * (T_surface - T★)
end

@inline function surface_salinity_relaxation(i, j, grid, clock, fields, p)
    time = clock.time

    n₁ = current_time_index(time, Nmonths)
    n₂ = next_time_index(time, Nmonths)

    @inbounds begin
        S★₁ = p.S★[i, j, n₁]
        S★₂ = p.S★[i, j, n₂]
        S_surface = fields.S[i, j, grid.Nz]
    end

    S★ = cyclic_interpolate(S★₁, S★₂, time)
                                
    return p.λ * (S_surface - S★)
end

T_surface_relaxation_bc = FluxBoundaryCondition(surface_temperature_relaxation,
                                                discrete_form = true,
                                                parameters = (λ = Δz_top/7days, T★ = target_sea_surface_temperature))

S_surface_relaxation_bc = FluxBoundaryCondition(surface_salinity_relaxation,
                                                discrete_form = true,
                                                parameters = (λ = Δz_top/7days, S★ = target_sea_surface_salinity))

u_bcs = FieldBoundaryConditions(top = u_wind_stress_bc, bottom = u_bottom_drag_bc)
v_bcs = FieldBoundaryConditions(top = v_wind_stress_bc, bottom = v_bottom_drag_bc)
T_bcs = FieldBoundaryConditions(top = T_surface_relaxation_bc)
S_bcs = FieldBoundaryConditions(top = S_surface_relaxation_bc)

free_surface = ImplicitFreeSurface(solver_method=:HeptadiagonalIterativeSolver, preconditioner_method = :SparseInverse,
                                   preconditioner_settings = (ε = 0.01, nzrel = 10))

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState())

model = HydrostaticFreeSurfaceModel(grid = grid,
                                    free_surface = free_surface,
				                    momentum_advection = VectorInvariant(),  
                                    coriolis = HydrostaticSphericalCoriolis(scheme = EnergyConservingScheme()),
                                    buoyancy = buoyancy,
                                    tracers = (:T, :S),
                                    closure = (horizontal_diffusivity, vertical_diffusivity, convective_adjustment, biharmonic_viscosity),
                                    boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs),
                                    forcing = (u=Fu, v=Fv),
                                    tracer_advection = WENO5(grid = underlying_grid))

#####
##### Initial condition:
#####

u, v, w = model.velocities
η = model.free_surface.η
T = model.tracers.T
S = model.tracers.S

file_init = jldopen("initial_conditions-1degree.jld2")

@info "Reading initial conditions"
T_init = file_init["T"]
S_init = file_init["S"]

set!(model, T=T_init, S=S_init)
fill_halo_regions!(T)
fill_halo_regions!(S)

@info "model initialized"

#####
##### Simulation setup
#####

Δt = 6minutes 

simulation = Simulation(model, Δt = Δt, stop_time = Nyears*years)

start_time = [time_ns()]

function progress(sim)
    wall_time = (time_ns() - start_time[1]) * 1e-9

    η = model.free_surface.η
    u = model.velocities.u
    @info @sprintf("Time: % 12s, iteration: %d, max(|η|): %.2e m, max(|u|): %.2e ms⁻¹, wall time: %s",
                    prettytime(sim.model.clock.time),
                    sim.model.clock.iteration,
                    maximum(abs, η), maximum(abs, u),
                    prettytime(wall_time))

    start_time[1] = time_ns()

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))


#=
u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
η = model.free_surface.η

output_fields = (; u, v, T, S, η)
save_interval = 5days

u2 = Field(u * u)
v2 = Field(v * v)
w2 = Field(w * w)
η2 = Field(η * η)
T2 = Field(T * T)

outputs = (; u, v, T, S, η)
average_outputs = (; u, v, T, S, η, u2, v2, T2, η2)

simulation.output_writers[:surface_fields] = JLD2OutputWriter(model, (; u, v, T, S, η),
                                                              schedule = TimeInterval(save_interval),
                                                              prefix = output_prefix * "_surface",
                                                              indices = (:, :, grid.Nz), 
                                                              force = true)

simulation.output_writers[:averages] = JLD2OutputWriter(model, average_outputs,
                                                              schedule = AveragedTimeInterval(4*30days, window=4*30days),
                                                              prefix = output_prefix * "_averages",
                                                              force = true)


simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                        schedule = TimeInterval(6*30days),
                                                        prefix = output_prefix * "_checkpoint",
                                                        force = true)
=#
# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

run!(simulation, pickup = pickup_file)

@info """

    Simulation took $(prettytime(simulation.run_wall_time))
    Background diffusivity: $background_diffusivity
    Free surface: $(typeof(model.free_surface).name.wrapper)
    Time step: $(prettytime(Δt))
"""

