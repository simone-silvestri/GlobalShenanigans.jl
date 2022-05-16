using GLMakie 

heatmap(randn(2,2))
contour(randn(2,2))
lines(randn(2))

fig = Figure()
ax = Axis(fig[1,1])
hm = heatmap!(ax, randn(2,2))
Colorbar(fig[2, 2], hm)