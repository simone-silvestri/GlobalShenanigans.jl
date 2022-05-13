using GLMakie, Printf

uvel = Array(interior(simulation.model.velocities.u))
vvel = Array(interior(simulation.model.velocities.v))
Ttrac = Array(interior(simulation.model.tracers.T))
Strac = Array(interior(simulation.model.tracers.S))

hc = buoyancy.equation_of_state.haline_contraction
ac = buoyancy.equation_of_state.thermal_expansion
ρ_prime = hc * Strac - ac * Ttrac

z = Array(grid.grid.zᵃᵃᶜ[1:end-4])
λc = Array(grid.grid.λᶜᵃᵃ[1:end-4])
# λf = Array(grid.grid.λᶠᵃᵃ[3:end-4])  # never used since periodic
φc = Array(grid.grid.φᵃᶜᵃ[1:end-4])
φf = Array(grid.grid.φᵃᶠᵃ[1:end-4])

continentsu = (uvel .== 0.0)
continentsv = (vvel .== 0.0)
continentst = (Ttrac .== 0.0)
uvel[continentsu] .= NaN
vvel[continentsv] .= NaN
Ttrac[continentst] .= NaN
Strac[continentst] .= NaN
ρ_prime[continentst] .= NaN

fig = Figure(resolution=(1400, 1200))
options = (; titlesize=30)
ustring = "U [m/s]"
vstring = "V [m/s]"
bstring = "T [Cᵒ]"
tstring = "S [psu]"

ax1 = Axis(fig[2, 1]; title=ustring, options...)
ax2 = Axis(fig[2, 3]; title=vstring, options...)
ax3 = Axis(fig[3, 1]; title=bstring, options...)
ax4 = Axis(fig[3, 3]; title=tstring, options...)

maxind = size(uvel)[end]
surface_slider = Slider(fig[4, 1:4], range=1:1:maxind, startvalue=maxind)
si = surface_slider.value # latitude index, lati for short

zstring = @lift(@sprintf("%0.1e", z[$si]))

ax0 = Label(fig[1, 1:4], text= @lift("Fields at z=" * $zstring * "[m]"), textsize=40)
hm1 = heatmap!(ax1, λc, φc, @lift(uvel[:, :, $si]), colormap=:balance, colorrange=(-0.2, 0.2), nan_color=:black, interpolate=true)
hm2 = heatmap!(ax2, λc, φf, @lift(vvel[:, :, $si]), colormap=:balance, colorrange=(-0.2, 0.2), nan_color=:black, interpolate=true)
hm3 = heatmap!(ax3, λc, φc, @lift(Ttrac[:, :, $si]), colormap=:thermometer, nan_color=:black, interpolate=true)
hm4 = heatmap!(ax4, λc, φc, @lift(Strac[:, :, $si]), colormap=:thermometer, nan_color=:black, interpolate=true)

Colorbar(fig[2, 2], hm1, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)
Colorbar(fig[3, 2], hm3, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)
Colorbar(fig[2, 4], hm2, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)
Colorbar(fig[3, 4], hm4, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)

for ax in [ax1, ax2, ax3, ax4]
    ax.limits = (extrema(λc)..., extrema(φc)...)

    ax.xlabel = "Longitude [ᵒ]"
    ax.ylabel = "Latitude [ᵒ]"
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.yticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])
    ax.xticks = ([-160, -120, -80, -40, 0, 40, 80, 120, 160], ["160W", "120W", "80W", "40W", "0", "40E", "80E", "120E", "160E"])
end

##viz
fig2 = Figure(resolution=(1400, 1200))
ax5 = Axis(fig2[2, 1]; title=ustring, options...)
ax6 = Axis(fig2[2, 3]; title=vstring, options...)
ax7 = Axis(fig2[3, 1]; title=bstring, options...)
ax8 = Axis(fig2[3, 3]; title=tstring, options...)

maxind = size(Ttrac)[2]
lat_slider = Slider(fig2[4, 1:4], range=1:1:maxind, startvalue=1)
lati = lat_slider.value # latitude index, lati for short

plati = @lift($lati < 75 ? string(75 - $lati + 1) * "S" : string($lati - 75) * "N")
latstring = plati
ax01 = Label(fig2[1, 1:4], text=@lift("Fields at latitude=" * $latstring), textsize=40)

hm5 = heatmap!(ax5, λc, z, @lift(uvel[:, $lati, :]), colormap=:balance, colorrange=(-0.1, 0.1), nan_color=:black, interpolate=true)
hm6 = heatmap!(ax6, λc, z, @lift(vvel[:, $lati, :]), colormap=:balance, colorrange=(-0.05, 0.05), nan_color=:black, interpolate=true)
hm7 = heatmap!(ax7, λc, z, @lift(Ttrac[:, $lati, :]), colormap=:thermometer, nan_color=:black, interpolate=true)
hm8 = heatmap!(ax8, λc, z, @lift(Strac[:, $lati, :]), colormap=:thermometer, nan_color=:black, interpolate=true)

Colorbar(fig2[2, 2], hm5, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)
Colorbar(fig2[3, 2], hm7, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)
Colorbar(fig2[2, 4], hm6, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)
Colorbar(fig2[3, 4], hm8, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)

for ax in [ax5, ax6, ax7, ax8]
    contour!(ax, λc, z, @lift(ρ_prime[:, $lati, :]), nan_color=:black, color=:black, levels=20)
end

for ax in [ax5, ax6, ax7, ax8]
    ax.limits = (extrema(λc)..., extrema(z)...)

    ax.xlabel = "Longitude [ᵒ]"
    ax.ylabel = "Depth [km]"
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.yticks = ([0, -1000, -2000, -3000, -4000, -5000], ["0", "1", "2", "3", "4", "5"])
    ax.xticks = ([-160, -120, -80, -40, 0, 40, 80, 120, 160], ["160W", "120W", "80W", "40W", "0", "40E", "80E", "120E", "160E"])
end

## viz 
fig3 = Figure(resolution=(1400, 1200))
ax9 = Axis(fig3[2, 1]; title=ustring, options...)
ax10 = Axis(fig3[2, 3]; title=vstring, options...)
ax11 = Axis(fig3[3, 1]; title=bstring, options...)
ax12 = Axis(fig3[3, 3]; title=tstring, options...)

maxind = size(Ttrac)[1]
long_slider = Slider(fig3[4, 1:4], range=1:1:maxind, startvalue=1)
longi = long_slider.value # latitude index, lati for short
# println("looking at longitude ", longi)
longstring = @lift($longi < 180 ? string(180 - $longi) * "W" : string($longi - 180 ) * "E")
ax01 = Label(fig3[1, 1:4], text=@lift("Fields at longitude=" * $longstring), textsize=40)
hm9 = heatmap!(ax9, φc, z, @lift(uvel[$longi, :, :]), colormap=:balance, colorrange=(-0.1, 0.1), nan_color=:black, interpolate=true)
hm10 = heatmap!(ax10, φf, z, @lift(vvel[$longi, :, :]), colormap=:balance, colorrange=(-0.05, 0.05), nan_color=:black, interpolate=true)
hm11 = heatmap!(ax11, φc, z, @lift(Ttrac[$longi, :, :]), colormap=:thermometer, nan_color=:black, interpolate=true)
hm12 = heatmap!(ax12, φc, z, @lift(Strac[$longi, :, :]), colormap=:thermometer, nan_color=:black, interpolate=true)

Colorbar(fig3[2, 2], hm9, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)
Colorbar(fig3[3, 2], hm11, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)
Colorbar(fig3[2, 4], hm10, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)
Colorbar(fig3[3, 4], hm12, height=Relative(3 / 4), width=25, ticklabelsize=30,
    labelsize=30, ticksize=25, tickalign=1,)

for ax in [ax9, ax10, ax11, ax12]
    contour!(ax, φc, z, @lift(ρ_prime[$longi, :, :]), nan_color=:black, color=:black, levels=20)
end


for ax in [ax9, ax10, ax11, ax12]
    ax.limits = (extrema(φc)..., extrema(z)...)

    ax.xlabel = "Latitude [ᵒ]"
    ax.ylabel = "Depth [km]"
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.yticks = ([0, -1000, -2000, -3000, -4000, -5000], ["0", "1", "2", "3", "4", "5"])
    ax.xticks = ([-80, -60, -30, 0, 30, 60, 80], ["80S", "60S", "30S", "0", "30N", "60N", "80N"])
end

##
iterator = collect(long_slider.attributes.range[][1]:10:long_slider.attributes.range[][end])
framerate = 10
record(fig3, "fixed_longitude.mp4", iterator;
    framerate=framerate) do its
    longi[] = its
end