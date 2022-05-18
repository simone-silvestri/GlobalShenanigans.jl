using GLMakie


fig = Figure()
ax = Axis(fig[1, 1])
ind1 = 170
ind2 = 89
# upper_quantile, median_quantile, lower_quantile = return_statistics(ecco_temperature)
xx = median_quantile[170, 89, :]
yy = ecco_grid.z # -log.(-ecco_grid.z)
lines!(ax, xx, yy, color=:black)
# errorbars!(ax, xx, yy, xx - lower_quantile[ind1, ind2, :], upper_quantile[ind1, ind2, :] - xx, direction=:x, whiskerwidth=3)
# band!(ax, xx, yy, xx - lower_quantile[ind1, ind2, :], upper_quantile[ind1, ind2, :] - xx, color = (:blue, 0.2))
l1 = [(upper_quantile[ind1, ind2, i], yy[i]) for i in 1:50]
l2 = [(lower_quantile[ind1, ind2, i], yy[i]) for i in 1:50]
poly!(ax, Point2f[l1..., reverse(l2)...], color=(:blue, 0.1))

ylims!(ax, s(-100, 0))
