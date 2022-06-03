using Statistics
using JLD2
using Printf
using GLMakie

using Dates: AbstractTime

# 0.25 degree resolution
function heatmap_surface_video_variable(folder, size)

    Nx, Ny, Nz = size

    output_prefix = "near_global_lat_lon_$(Nx)_$(Ny)_$(Nz)_fine"
    surface_file = jldopen(folder * output_prefix * "_surface.jld2")
    vorticity    = jldopen(folder * "vorticity.jld2")

    iterations = parse.(Int, keys(surface_file["timeseries/t"]))

    iterations = iterations[1:end-1]
    iter = Observable(iterations[1])

    ki(iter) = surface_file["timeseries/v/" * string(iter)][:, 1:end-1, 1].^2 .+ surface_file["timeseries/u/" * string(iter)][:, :, 1].^2
    ti(iter) = string(surface_file["timeseries/t/" * string(iter)] / 3600 /24)

    diter = iterations[3] - iterations[2]
 
    ζi(iter) = vorticity["vort"][:, 2:end, Int(floor(iter/diter)) + 1]

    bathy = surface_file["timeseries/η/" * string(iterations[end-1])][:, :, 1]
    bathyv = surface_file["timeseries/v/" * string(iterations[end-1])][:, :, 1]

    bathy[bathy .== 0]    .= NaN
    bathy[(!).(isnan.(bathy))] .= 0
    bathyv[bathyv .== 0]    .= NaN
    bathyv[(!).(isnan.(bathyv))] .= 0

    k = @lift ki($iter) .+ bathy
    ζ = @lift ζi($iter) .+ bathy

    max_η = 1.5
    min_η = - max_η
    max_u = 0.3
    min_u = - max_u

    fig = Figure(resolution = (3500, 1000))
    ga  = fig[1, 1] = GridLayout() 

    ax = Axis(ga[1, 1]) #, title="Surface vertical vorticity (s⁻¹)")
    hm = GLMakie.heatmap!(ax, ζ, colorrange=(-1e-5, 1e-5), colormap=:blues, nan_color=:black, interpolate = true)
    hidedecorations!(ax, grid = false)

    ax = Axis(ga[1, 2]) #, title="Surface kinetic energy (m²s⁻²)")
    hm = GLMakie.heatmap!(ax, k, colorrange=(0.0, 2e-1), colormap=:solar, nan_color=:black, interpolate = true)
    hidedecorations!(ax, grid = false)
    colgap!(ga, 0)
    # cb = Colorbar(fig[1, 4], hm)

    # title_str = @lift "Earth day = " * ti($iter)
    # ax_t = fig[0, :] = Label(fig, title_str)

    GLMakie.record(fig, output_prefix * ".mp4", iterations[1:end-3], framerate=15) do i
        @info "Plotting iteration $i of $(iterations[end])..."
        iter[] = i
    end

    display(fig)
end


function heatmap_surface_video_variable(folder1, size1, folder2, size2)

    Nx1, Ny1, Nz1 = size1
    Nx2, Ny2, Nz2 = size2

    output_prefix1 = "near_global_lat_lon_$(Nx1)_$(Ny1)_$(Nz1)_fine"
    surface_file1  = jldopen(folder1 * output_prefix1 * "_surface.jld2")
    vorticity1     = jldopen(folder1 * "vorticity.jld2")

    output_prefix2 = "near_global_lat_lon_$(Nx2)_$(Ny2)_$(Nz2)_fine"
    surface_file2  = jldopen(folder2 * output_prefix2 * "_surface.jld2")
    vorticity2     = jldopen(folder2 * "vorticity.jld2")

    iterations1 = parse.(Int, keys(surface_file1["timeseries/t"]))
    iterations2 = parse.(Int, keys(surface_file2["timeseries/t"]))

    if length(iterations2) > length(iterations1)
        iterations = iterations1
    else
        iterations = iterations2
    end

    iter = Observable(1)

    ki1(iter) = surface_file1["timeseries/v/" * string(iterations1[iter])][:, 1:end-1, 1].^2 .+ surface_file1["timeseries/u/" * string(iterations1[iter])][:, :, 1].^2
    ki2(iter) = surface_file2["timeseries/v/" * string(iterations2[iter])][:, 1:end-1, 1].^2 .+ surface_file2["timeseries/u/" * string(iterations2[iter])][:, :, 1].^2
 
    ζi1(iter) = vorticity1["vort"][:, 2:end, iter]
    ζi2(iter) = vorticity2["vort"][:, 2:end, iter]

    bathy1 = vorticity1["vort"][:, 2:end, end]
    bathy2 = vorticity2["vort"][:, 2:end, end]

    bathy1[bathy1 .== 0]    .= NaN
    bathy1[(!).(isnan.(bathy1))] .= 0
    bathy2[bathy2 .== 0]    .= NaN
    bathy2[(!).(isnan.(bathy2))] .= 0

    k1 = @lift ki1($iter) .+ bathy1
    ζ1 = @lift ζi1($iter) .+ bathy1
    k2 = @lift ki2($iter) .+ bathy2
    ζ2 = @lift ζi2($iter) .+ bathy2

    fig = Figure(resolution = (3500, 2000))
    ga  = fig[1, 1] = GridLayout() 

    ax = Axis(ga[1, 1]) 
    hm = GLMakie.heatmap!(ax, ζ1, colorrange=(-1e-5, 1e-5), colormap=:blues, nan_color=:black, interpolate = true)
    hidedecorations!(ax, grid = false)

    ax = Axis(ga[1, 2]) 
    hm = GLMakie.heatmap!(ax, k1, colorrange=(0.0, 2e-1), colormap=:solar, nan_color=:black, interpolate = true)
    hidedecorations!(ax, grid = false)

    ax = Axis(ga[2, 1]) 
    hm = GLMakie.heatmap!(ax, ζ2, colorrange=(-1e-5, 1e-5), colormap=:blues, nan_color=:black, interpolate = true)
    hidedecorations!(ax, grid = false)

    ax = Axis(ga[2, 2]) 
    hm = GLMakie.heatmap!(ax, k2, colorrange=(0.0, 2e-1), colormap=:solar, nan_color=:black, interpolate = true)
    hidedecorations!(ax, grid = false)
    
    colgap!(ga, 0)
    rowgap!(ga, 0)

    GLMakie.record(fig, "all-video.mp4", 1:length(iterations), framerate=15) do i
        @info "Plotting iteration $i of $(length(iterations))..."
        iter[] = i
    end

    display(fig)
end


using GMT

function gmt_image(ki, lim, dlim, λ, ϕ, colormap, figname)
    k  = ki(iterations[end]) .+ bathy
    k2 = Array(k').^(0.5)
    k2 = clamp.(k2, lim[1], lim[2])
    G2 = mat2grid(k2, x = λ, y = ϕ) 
    topo = makecpt(range=(lim[1], lim[2]+dlim, 0.001); continuous=true, colormap)
    grdimage(G2; show = true, projection = :Robinson, region = (-180, 180, -75, 75), color = topo, figname)
end


k  = ki(iterations[end]) .+ bathy
k2 = Array(k')
k3 = clamp.(k2, 0.0, 0.47)
G2 = mat2grid(k3, x = λ, y = ϕ) 
topo = makecpt(range=(0.0, 0.5, 0.001); continuous=true)
grdimage(G2, show = true, projection = :Robinson, region = (-180, 180, -75, 75), color = topo, figname="test2.png")

gmt_image(ki, (0.0, 0.5), 0.001, λ, ϕ, :solar, "test.png")