using GLMakie

function closest_index(value, field)
    argmin(abs.(value .- field))
end

function minmax_extrema(coord1, coord2)
    maxmincoord = max(minimum(coord1), minimum(coord2))
    minmaxcoord = min(maximum(coord1), maximum(coord2))
    return (maxmincoord, minmaxcoord)
end

function plot_extrema(grid1, grid2)
    λ1, ϕ1, z1 = grid1
    λ2, ϕ2, z2 = grid2

    λ = minmax_extrema(λ1, λ2)
    ϕ = minmax_extrema(ϕ1, ϕ2)
    z = minmax_extrema(z1, z2)

    return (; λ, ϕ, z)
end

function comparison(grid1, solution1, grid2, solution2; quantiles=(0.01, 0.99), colormap=:linear_protanopic_deuteranopic_kbw_5_98_c40_n256, quick=nothing)
    if isnothing(quick)
        nothing
    elseif quick == :linear
        colormap = :linear_protanopic_deuteranopic_kbw_5_98_c40_n256
    elseif quick == :diverging
        colormap = :diverging_protanopic_deuteranopic_bwy_60_95_c32_n256
    else
        println("quickcolormap option does not exist")
    end

    fig = Figure()
    ax1 = Axis(fig[2, 2])
    ax2 = Axis(fig[2, 3])
    ax3 = Axis(fig[3, 2:3])

    λ1, ϕ1, z1 = grid1
    λ2, ϕ2, z2 = grid2

    p_extrema = plot_extrema(grid1, grid2)

    λ_slider = Slider(fig[4, 2:3], range=-180:180, startvalue=0)
    λvalue = λ_slider.value
    λ1_index = @lift(closest_index($λvalue, λ1))
    λ2_index = @lift(closest_index($λvalue, λ2))

    ϕ_slider = Slider(fig[3, 1], range=-75:75, startvalue=0, horizontal=false)
    ϕvalue = ϕ_slider.value
    ϕ1_index = @lift(closest_index($ϕvalue, ϕ1))
    ϕ2_index = @lift(closest_index($ϕvalue, ϕ2))

    field1 = @lift(solution1[$λ1_index, :, :])
    field2 = @lift(solution2[$λ2_index, :, :])

    surface_index = closest_index(0.0, z1) # 0 is assumed to be the surface
    surface_field = solution1[:, :, surface_index]

    colorrange = quantile.(Ref(filter(!isnan, solution1)), quantiles)
    surfacecolorrange = quantile.(Ref(filter(!isnan, surface_field)), quantiles)

    hm = heatmap!(ax1, ϕ1, z1, field1, colorrange=colorrange, colormap=colormap)
    heatmap!(ax2, ϕ2, z2, field2, colorrange=colorrange, colormap=colormap)
    heatmap!(ax3, λ1, ϕ1, surface_field, colorrange=surfacecolorrange, colormap=colormap)
    vl = vlines!(ax3, @lift(λ1[$λ1_index]), color=:orange, linewidth=3)
    hl = hlines!(ax3, @lift(ϕ1[$ϕ1_index]), color=:yellow, linewidth=3)

    Colorbar(fig[2, 4], hm, height=Relative(3 / 4), width=25, ticklabelsize=30,
        labelsize=30, ticksize=25, tickalign=1,)

    for ax in [ax1, ax2]
        ax.limits = (p_extrema.ϕ..., p_extrema.z...)

        ax.xlabel = "Latitude [ᵒ]"
        ax.ylabel = "Depth [km]"
        ax.xlabelsize = 25
        ax.ylabelsize = 25

        ax.xticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])
        ax.yticks = ([0, -1000, -2000, -3000, -4000, -5000], ["0", "1", "2", "3", "4", "5"])
    end

    ax3.limits = (p_extrema.λ..., p_extrema.ϕ...,)
    ax3.xlabel = "Longitude [ᵒ]"
    ax3.ylabel = "Latitude [ᵒ]"
    ax3.xlabelsize = 25
    ax3.ylabelsize = 25
    ax3.xticks = ([-160, -120, -80, -40, 0, 40, 80, 120, 160], ["160W", "120W", "80W", "40W", "0", "40E", "80E", "120E", "160E"])
    ax3.yticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])

    return fig
end