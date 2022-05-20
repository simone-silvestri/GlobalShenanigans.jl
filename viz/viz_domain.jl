using Statistics
using JLD2
using Printf
using GLMakie

using Dates: AbstractTime

# 0.25 degree resolution
Nx = 1440
Ny = 600
Nz = 48

folder = "./"

output_prefix = folder * "near_global_lat_lon_$(Nx)_$(Ny)_$(Nz)_fine"

surface_file = jldopen(output_prefix * "_surface.jld2")
vorticity    = jldopen(folder * "vorticity.jld2")

iterations = parse.(Int, keys(surface_file["timeseries/t"]))

iterations = iterations[1:end-1]
iter = Observable(iterations[1])

ηi(iter) = surface_file["timeseries/η/" * string(iter)][:, :, 1]
ui(iter) = - surface_file["timeseries/u/" * string(iter)][:, :, 1]
vi(iter) = surface_file["timeseries/v/" * string(iter)][:, :, 1]
ki(iter) = surface_file["timeseries/v/" * string(iter)][:, 1:end-1, 1].^2 .+ surface_file["timeseries/u/" * string(iter)][:, :, 1].^2
Ti(iter) = surface_file["timeseries/T/" * string(iter)][:, :, 1]
ti(iter) = string(surface_file["timeseries/t/" * string(iter)] / 3600 /24)

diter = iterations[3] - iterations[2]
 
ζi(iter) = vorticity["vort"][:, 2:end, Int(floor(iter/diter)) + 1]

bathy = surface_file["timeseries/η/" * string(iterations[end-1])][:, :, 1]
bathyv = surface_file["timeseries/v/" * string(iterations[end-1])][:, :, 1]

bathy[bathy .== 0]    .= NaN
bathy[(!).(isnan.(bathy))] .= 0
bathyv[bathyv .== 0]    .= NaN
bathyv[(!).(isnan.(bathyv))] .= 0

η = @lift ηi($iter) .+ bathy
u = @lift ui($iter) .+ bathy
k = @lift ki($iter) .+ bathy
v = @lift vi($iter) .+ bathyv
T = @lift Ti($iter) .+ bathy
ζ = @lift ζi($iter) .+ bathy

max_η = 1.5
min_η = - max_η
max_u = 0.3
min_u = - max_u
max_T = 32
min_T = 0

# fig = Figure(resolution = (3500, 1000))
# 
# ax = Axis(fig[1, 1], title="Free surface displacement (m)")
# colormap = :hot; colorrange = (-2, 1); contour_levels = collect(-0.7:0.1:0.9)
# hm    = GLMakie.heatmap!(ax, η, nan_color = :black, colorrange = colorrange, colormap = colormap, interpolate = true)
# cplot = GLMakie.contour!(ax, η, levels = contour_levels, color = :black, interpolate = true)
# cb = Colorbar(fig[1, 2], hm)
# 
# ax = Axis(fig[1, 3], title="Surface temperature (ᴼC)")
# contour_levels = collect(0:2:32)
# hm    = GLMakie.heatmap!(ax, T, colorrange=(min_T, max_T), colormap=:thermal) #, background_color = :black)
# cplot = GLMakie.contour!(ax, T, levels = contour_levels, color = :black, interpolate = true)
# cb    = Colorbar(fig[1, 4], hm)
# 
# ax = Axis(fig[2, 1], title="East-west surface velocity (m s⁻¹)")
# hm = GLMakie.heatmap!(ax, u, colorrange=(min_u, max_u), colormap=:balance, nan_color = :black) #, background_color = :black)
# cb = Colorbar(fig[2, 2], hm)
# 
# ax = Axis(fig[1, 3], title="North-south surface velocity (m s⁻¹)")
# hm = GLMakie.heatmap!(ax, v, colorrange=(min_u, max_u), colormap=:balance, nan_color = :black) #, background_color = :black)
# # cb = Colorbar(fig[1, 4], hm)
# 
# ax = Axis(fig[1, 1], title="Surface vertical vorticity (s⁻¹)")
# hm = GLMakie.heatmap!(ax, ζ, colorrange=(-1e-5, 1e-5), colormap=:blues, nan_color=:black, interpolate = true)
# cb = Colorbar(fig[1, 2], hm)
# 
# ax = Axis(fig[1, 3], title="Surface kinetic energy (m²s⁻²)")
# hm = GLMakie.heatmap!(ax, k, colorrange=(0.0, 2e-1), colormap=:solar, nan_color=:black, interpolate = true)
# cb = Colorbar(fig[1, 4], hm)
# 
# # ax = Axis(fig[1, 1], title="North-south surface velocity (m s⁻¹) @ -900 m")
# # hm = GLMakie.heatmap!(ax, u, colorrange=(min_u, max_u), colormap=:balance, nan_color=:black)
# # cb = Colorbar(fig[1, 2], hm)
# 
# title_str = @lift "Earth day = " * ti($iter)
# ax_t = fig[0, :] = Label(fig, title_str)
# 
# GLMakie.record(fig, output_prefix * ".mp4", iterations[1:end-3], framerate=15) do i
#     @info "Plotting iteration $i of $(iterations[end])..."
#     iter[] = i
# end
# 
# display(fig)

