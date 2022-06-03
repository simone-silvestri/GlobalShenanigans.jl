using GLMakie, Statistics

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

function show_lon(λ)
    if λ < 0
        return "$(-λ)ᵒW"
    else
        return "$(λ)ᵒE"
    end
end

function show_lat(ϕ)
    if ϕ < 0
        return "$(-ϕ)ᵒS"
    else
        return "$(ϕ)ᵒN"
    end
end

function comparison(grid1, solution11, solution12, grid2, solution21, solution22, grid3, solution31, solution32; quantiles=(0.01, 0.99), colormap=:linear_protanopic_deuteranopic_kbw_5_98_c40_n256, quick=nothing)
    if isnothing(quick)
        nothing
    elseif quick == :linear
        colormap = :linear_protanopic_deuteranopic_kbw_5_98_c40_n256
    elseif quick == :diverging
        colormap = :diverging_protanopic_deuteranopic_bwy_60_95_c32_n256
    else
        println("quickcolormap option does not exist")
    end

    fig = Figure(resolution = (1664, 1053))
    ax1 = Axis(fig[2, 2])
    ax2 = Axis(fig[2, 3])
    ax3 = Axis(fig[3, 2:3])

    λ1, ϕ1, z1 = grid1
    λ2, ϕ2, z2 = grid2
    λ3, ϕ3, z3 = grid3

    p_extrema = plot_extrema(grid1, grid2)

    λ_slider = Slider(fig[4, 2:3], range=-180:180, startvalue=0)
    λvalue = λ_slider.value
    λ1_index = @lift(closest_index($λvalue, λ1))
    λ2_index = @lift(closest_index($λvalue, λ2))
    λ3_index = @lift(closest_index($λvalue, λ3))

    ϕ_slider = Slider(fig[3, 1], range=-75:75, startvalue=0, horizontal=false)
    ϕvalue = ϕ_slider.value
    ϕ1_index = @lift(closest_index($ϕvalue, ϕ1))
    ϕ2_index = @lift(closest_index($ϕvalue, ϕ2))
    ϕ3_index = @lift(closest_index($ϕvalue, ϕ3))

    ϕstring = @lift(show_lat($ϕvalue))
    λstring = @lift(show_lon($λvalue))
    ax4 = Axis(fig[2, 5], xlabel = @lift("Temperature [ᵒC] at Longitude=" * $λstring * " and Latitude=" * $ϕstring))
    ax5 = Axis(fig[3, 5], xlabel = @lift("Salinity [psu] at Longitude=" * $λstring * " and Latitude=" * $ϕstring))
    ax0 = Label(fig[1, 1:4], text=@lift("Temperature [ᵒC] at Longitude=" * $λstring * ",    0.25ᵒ setup (left) and 1ᵒ setup (right)"), textsize=20)

    field1 = @lift(solution11[$λ1_index, :, :])
    field2 = @lift(solution21[$λ2_index, :, :])

    field11_profile = @lift(solution11[$λ1_index, $ϕ1_index, :])
    field21_profile = @lift(solution21[$λ2_index, $ϕ2_index, :])
    field31_profile = @lift(solution31[$λ3_index, $ϕ3_index, :])    
    
    field12_profile = @lift(solution12[$λ1_index, $ϕ1_index, :])
    field22_profile = @lift(solution22[$λ2_index, $ϕ2_index, :])
    field32_profile = @lift(solution32[$λ3_index, $ϕ3_index, :])

    
    surface_index = closest_index(0.0, z1) # 0 is assumed to be the surface
    surface_field = solution11[:, :, surface_index]

    colorrange1 = quantile.(Ref(filter(!isnan, solution11)), quantiles)
    colorrange2 = quantile.(Ref(filter(!isnan, solution12)), quantiles)
    surfacecolorrange = quantile.(Ref(filter(!isnan, surface_field)), quantiles)

    hm = heatmap!(ax1, ϕ1, z1, field1, colorrange=colorrange1, colormap=colormap, nan_color=:black)
    heatmap!(ax2, ϕ2, z2, field2, colorrange=colorrange1, colormap=colormap, nan_color=:black)
    heatmap!(ax3, λ1, ϕ1, surface_field, colorrange=surfacecolorrange, colormap=colormap, nan_color=:black)
    vl = vlines!(ax3, @lift(λ1[$λ1_index]), color=:orange, linewidth=3)
    hl = hlines!(ax3, @lift(ϕ1[$ϕ1_index]), color=:yellow, linewidth=3)

    line1 = lines!(ax4, field11_profile, z1, color = :blue , label = "0.25ᵒ simulation")
    line2 = lines!(ax4, field21_profile, z2, color = :red  , label = "1ᵒ simulation")
    line3 = lines!(ax4, field31_profile, z3, color = :green, label = "interpolated observations")
    axislegend(ax4, position = :rb)

    line1 = lines!(ax5, field12_profile, z1, color = :blue , label = "0.25ᵒ simulation")
    line2 = lines!(ax5, field22_profile, z2, color = :red  , label = "1ᵒ simulation")
    line3 = lines!(ax5, field32_profile, z3, color = :green, label = "interpolated observations")
    axislegend(ax5, position = :rb)

    Colorbar(fig[2, 4], hm, height=Relative(3 / 4), width=15, ticklabelsize=15,
        labelsize=30, ticksize=15, tickalign=1,)

    for ax in [ax1, ax2]
        ax.limits = (p_extrema.ϕ..., p_extrema.z...)

        ax.xlabel = "Latitude [ᵒ]"
        ax.ylabel = "Depth [km]"
        ax.xlabelsize = 15
        ax.ylabelsize = 15

        ax.xticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])
        ax.yticks = ([0, -1000, -2000, -3000, -4000, -5000], ["0", "1", "2", "3", "4", "5"])
    end

    ax3.limits = (p_extrema.λ..., p_extrema.ϕ...,)
    ax3.xlabel = "Longitude [ᵒ]"
    ax3.ylabel = "Latitude [ᵒ]"
    ax3.xlabelsize = 15
    ax3.ylabelsize = 15
    ax3.xticks = ([-160, -120, -80, -40, 0, 40, 80, 120, 160], ["160W", "120W", "80W", "40W", "0", "40E", "80E", "120E", "160E"])
    ax3.yticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])

    ax4.limits = (colorrange1..., -4000, 0)
    ax4.yticks = ([0, -1000, -2000, -3000, -4000, -5000], ["0", "1", "2", "3", "4", "5"])
    ax4.ylabel = "Depth [km]"
    ax4.xlabelsize = 15
    ax4.ylabelsize = 15

    ax5.limits = (colorrange2..., -4000, 0)
    ax5.yticks = ([0, -1000, -2000, -3000, -4000, -5000], ["0", "1", "2", "3", "4", "5"])
    ax5.ylabel = "Depth [km]"
    ax5.xlabelsize = 15
    ax5.ylabelsize = 15

    ax5.limits = (33.5, 37.5, -4000, 0)
    display(fig)

    return fig, λ_slider.value, ϕ_slider.value
end

fig, λ, ϕ = comparison(grid, Ttrac, Strac, grid2, Ttrac2, Strac2, grid_ecco, Tecco, Secco)

path1 = ((30, -170 + i) for i in 0:50)
path2 = ((30 - i, -120) for i in 1:90)
path3 = ((-60, -120 +i) for i in 1:90)
path4 = ((-60 + i, -30) for i in 1:105)
path5 = ((45, -30  + i) for i in 1:10)
path6 = ((45 - i, -20 ) for i in 1:90)
path7 = ((-45, -20 + i) for i in 1:80)
path8 = ((-45 + i, 60 ) for i in 1:45)
path9 = ((0   , 60 + i) for i in 1:30)
path10= ((0 - i , 90  ) for i in 1:50)
path11= ((-50 , 90 + i) for i in 1:89)
path12= ((-50 , -179 +i) for i in 1:9)
path13= ((-50 +i, -170) for i in 1:80)

long_path = (path1..., path2..., path3..., path4..., path5..., path6..., path7..., path8..., path9..., path10..., path11..., path12..., path13...)

GLMakie.record(fig, "comparison-video.mp4", 1:length(long_path), framerate=10) do i
    @info "Plotting iteration $i..."
    ϕ[] = long_path[i][1]
    λ[] = long_path[i][2]
end