using GLMakie

function closest_index(value, field)
    argmin(abs.(value .- field))
end

function comparison(grid1, solution1, grid2, solution2)
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 2])
    λ1, ϕ1, z1 = grid1
    λ2, ϕ2, z2 = grid2

    λ_slider = Slider(fig[2, 1:2], range=-180:180, startvalue=0)
    λvalue = λ_slider.value # lonitude index, loni for short
    λ1_index = @lift(closest_index($λvalue, λ1))
    λ2_index = @lift(closest_index($λvalue, λ2))

    ϕ_slider = Slider(fig[1:2, 3], range=-75:75, startvalue=0, horizontal = false)
    ϕvalue = ϕ_slider.value # lonitude index, loni for short
    ϕ1_index = @lift(closest_index($ϕvalue, ϕ1))
    ϕ2_index = @lift(closest_index($ϕvalue, ϕ2))


    field1 = @lift(solution1[$λ1_index, :, :])
    field2 = @lift(solution2[$λ2_index, :, :])
    heatmap!(ax1, ϕ1, z1, field1)
    heatmap!(ax2, ϕ2, z2, field2)

    return fig
end