import Pkg
Pkg.activate("/home/ssilvest/Oceananigans.jl")

using Statistics
using JLD2
using Printf
using Plots
using Oceananigans
using Oceananigans.Units

using Oceananigans.Operators
using Oceananigans.MultiRegion
using Oceananigans.MultiRegion: multi_region_object_from_array
using Oceananigans.Fields: interpolate, Field
using Oceananigans.Architectures: arch_array
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis
using Oceananigans.BoundaryConditions
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom, inactive_node, peripheral_node
using CUDA: @allowscalar, device!
using Oceananigans.Operators
using Oceananigans.Operators: Δzᵃᵃᶜ
using Oceananigans: prognostic_fields

include("horizontal_visc.jl")

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

arch = GPU();
reference_density = 1029;

latitude = (-75, 75);

# 1/12 degree resolution
Nx = 4320;
Ny = 1800;
Nz = 48;

const Nyears  = 1;
const Nmonths = 12;
const thirty_days = 30days;

output_prefix = "near_global_lat_lon_$(Nx)_$(Ny)_$(Nz)_fine"
pickup_file   = false #"near_global_lat_lon_4320_1800_48_fine_checkpoint_iteration5000.jld2" 

#####
##### Load forcing files and inital conditions from ECCO version 4
##### https://ecco.jpl.nasa.gov/drive/files
##### Bathymetry is interpolated from ETOPO1 https://www.ngdc.noaa.gov/mgg/global/
#####

using DataDeps

path = "https://github.com/CliMA/OceananigansArtifacts.jl/raw/ss/new_hydrostatic_data_after_cleared_bugs/quarter_degree_near_global_input_data/"

dh = DataDep("quarter_degree_near_global_lat_lon",
    "Forcing data for global latitude longitude simulation",
    path * "z_faces-50-levels.jld2"
)

DataDeps.register(dh)

datadep"quarter_degree_near_global_lat_lon"

datadep_path = @datadep_str "quarter_degree_near_global_lat_lon/z_faces-50-levels.jld2"
file_z_faces = jldopen(datadep_path)

bathy_path = "../data/bathymetry-ad-hoc.jld2" # "smooth-bathymetry.jld2" #
bathymetry = jldopen(bathy_path)["bathymetry"]

τˣ = zeros(Nx, Ny, Nmonths)
τʸ = zeros(Nx, Ny, Nmonths)
T★ = zeros(Nx, Ny, Nmonths)
S★ = zeros(Nx, Ny, Nmonths)

path_bc = "../data/boundary_conditions-1-12degree.jld2"

# Files contain 1 year (1992) of 12 monthly averages
τˣ = jldopen(path_bc)["τˣ"] ./ reference_density;
τʸ = jldopen(path_bc)["τʸ"] ./ reference_density;
T★ = jldopen(path_bc)["Tˢ"];
S★ = jldopen(path_bc)["Sˢ"];

T_bounds = extrema(T★)
S_bounds = extrema(S★)

# Remember the convention!! On the surface a negative flux increases a positive decreases
bathymetry = arch_array(arch, bathymetry);

# Stretched faces taken from ECCO Version 4 (50 levels in the vertical)
z_faces = file_z_faces["z_faces"][3:end];

# A spherical domain
@show underlying_grid = LatitudeLongitudeGrid(arch,
                                              size = (Nx, Ny, Nz),
                                              longitude = (-180, 180),
                                              latitude = latitude,
                                              halo = (4, 4, 4),
                                              z = z_faces,
                                              precompute_metrics = true)

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry));

underlying_mrg = MultiRegionGrid(underlying_grid, partition = XPartition(4), devices = 4);
mrg            = MultiRegionGrid(grid,            partition = XPartition(4), devices = 4);

τˣ = multi_region_object_from_array(- τˣ, mrg);
τʸ = multi_region_object_from_array(- τʸ, mrg);

target_sea_surface_temperature = T★ = multi_region_object_from_array(T★, mrg);
target_sea_surface_salinity    = S★ = multi_region_object_from_array(S★, mrg);

#####
##### Physics and model setup
#####

using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: HorizontalDivergenceFormulation, HorizontalFormulation

νz = 5e-3
κz = 1e-4

biharmonic_viscosity   = smagorinsky_viscosity(HorizontalDivergenceFormulation(), grid; Cₛₘ = 4.5)
vertical_diffusivity   = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=νz, κ=κz)
convective_adjustment  = ConvectiveAdjustmentVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), convective_κz = 1.0)

closures = (biharmonic_viscosity, vertical_diffusivity, convective_adjustment)

#####
##### Boundary conditions / time-dependent fluxes 
#####

@inline current_time_index(time, tot_months)     = mod(unsafe_trunc(Int32, time / thirty_days), tot_months) + 1
@inline next_time_index(time, tot_months)        = mod(unsafe_trunc(Int32, time / thirty_days) + 1, tot_months) + 1
@inline cyclic_interpolate(u₁::Number, u₂, time) = u₁ + mod(time / thirty_days, 1) * (u₂ - u₁)

Δz_top    = @allowscalar Δzᵃᵃᶜ(1, 1, grid.Nz, grid.underlying_grid)
Δz_bottom = @allowscalar Δzᵃᵃᶜ(1, 1, 1, grid.underlying_grid)

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

u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τˣ);
v_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τʸ);

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1]
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1]

# Linear bottom drag:
μ = 0.001 # ms⁻¹

@inline is_immersed_drag_u(i, j, k, grid) = Int(peripheral_node(Face(), Center(), Center(), i, j, k-1, grid) & !inactive_node(Face(), Center(), Center(), i, j, k, grid))
@inline is_immersed_drag_v(i, j, k, grid) = Int(peripheral_node(Center(), Face(), Center(), i, j, k-1, grid) & !inactive_node(Center(), Face(), Center(), i, j, k, grid))                                

# Keep a constant linear drag parameter independent on vertical level
@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * is_immersed_drag_u(i, j, k, grid) * fields.u[i, j, k] / Δzᵃᵃᶜ(i, j, k, grid)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * is_immersed_drag_v(i, j, k, grid) * fields.v[i, j, k] / Δzᵃᵃᶜ(i, j, k, grid)

Fu = Forcing(u_immersed_bottom_drag, discrete_form = true, parameters = μ);
Fv = Forcing(v_immersed_bottom_drag, discrete_form = true, parameters = μ);

u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ);
v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ);

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
                                                parameters = (λ = Δz_top/7days, T★ = target_sea_surface_temperature));

S_surface_relaxation_bc = FluxBoundaryCondition(surface_salinity_relaxation,
                                                discrete_form = true,
                                                parameters = (λ = Δz_top/7days, S★ = target_sea_surface_salinity));

u_bcs = FieldBoundaryConditions(top = u_wind_stress_bc, bottom = u_bottom_drag_bc);
v_bcs = FieldBoundaryConditions(top = v_wind_stress_bc, bottom = v_bottom_drag_bc);
T_bcs = FieldBoundaryConditions(top = T_surface_relaxation_bc);
S_bcs = FieldBoundaryConditions(top = S_surface_relaxation_bc);

free_surface = ImplicitFreeSurface(solver_method=:HeptadiagonalIterativeSolver);

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState())

using Oceananigans.Advection: VelocityStencil

T_advection = WENO5(underlying_mrg, bounds = T_bounds)
S_advection = WENO5(underlying_mrg, bounds = S_bounds)

model = HydrostaticFreeSurfaceModel(grid = mrg,
                                    free_surface = free_surface,
				    momentum_advection = WENO5(vector_invariant = VelocityStencil()),
                                    coriolis = HydrostaticSphericalCoriolis(),
                                    buoyancy = buoyancy,
                                    tracers = (:T, :S),
				    closure = closures,
                                    boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs),
                                    forcing = (u=Fu, v=Fv),
				    tracer_advection = (T = T_advection, S = S_advection))

#####
##### Initial condition:
#####

u, v, w = model.velocities
η = model.free_surface.η
T = model.tracers.T
S = model.tracers.S

@info "Reading initial conditions"
file_init = jldopen("../data/evolved-initial-conditions-120days.jld2")

@info "initializing model"
T_init = multi_region_object_from_array(file_init["T"], mrg);
S_init = multi_region_object_from_array(file_init["S"], mrg);
u_init = multi_region_object_from_array(file_init["u"], mrg);
v_init = multi_region_object_from_array(file_init["v"], mrg);
η_init = multi_region_object_from_array(file_init["η"], mrg);
set!(model, T=T_init, S=S_init, u=u_init, v=v_init, η=η_init)

@info "model initialized"

#####
##### Simulation setup
#####

Δt = 120  # for initialization, then we can go up to 6 minutes?

simulation = Simulation(model, Δt = Δt, stop_iteration = 40000)

start_time = [time_ns()]

using Oceananigans.Utils 
using Oceananigans.MultiRegion: reconstruct_global_field

function progress(sim)
    wall_time = (time_ns() - start_time[1]) * 1e-9

    u = sim.model.velocities.u

    @info @sprintf("Time: % 12s, iteration: %d, max(|u|): %.2e ms⁻¹, wall time: %s", 
                    prettytime(sim.model.clock.time),
                    sim.model.clock.iteration, maximum(abs, u),
                    prettytime(wall_time))

    start_time[1] = time_ns()

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                        schedule = IterationInterval(5000),
                                                        prefix = output_prefix * "_checkpoint",
                                                        overwrite_existing = true)

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
η = model.free_surface.η

output_fields = (; u, v, T, S, η)
save_interval = 5days

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

run!(simulation, pickup = pickup_file)

@info """

    Simulation took $(prettytime(simulation.run_wall_time))
    Free surface: $(typeof(model.free_surface).name.wrapper)
    Time step: $(prettytime(Δt))
"""


