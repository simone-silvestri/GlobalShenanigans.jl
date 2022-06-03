using Oceananigans.TurbulenceClosures
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, ℑxyᶜᶜᵃ, ζ₃ᶠᶠᶜ, div_xyᶜᶜᶜ

@inline Dₜ(i, j, k, grid, u, v) = ∂xᶠᶠᶜ(i, j, k, grid, v) + ∂yᶠᶠᶜ(i, j, k, grid, u)
@inline Dₛ(i, j, k, grid, u, v) = ∂xᶜᶜᶜ(i, j, k, grid, u) - ∂yᶜᶜᶜ(i, j, k, grid, v)
@inline Δᶜᶜᶜ(i, j, k, grid)    = min(Δxᶜᶜᶜ(i, j, k, grid), Δyᶜᶜᶜ(i, j, k, grid))

function smagorinsky_viscosity(formulation, grid; Cₛₘ = 4.0)

    dx_min    = min_Δx(grid.underlying_grid)
    dy_min    = min_Δy(grid.underlying_grid)
    timescale = 10days

    @show C₄    = (Cₛₘ / π)^2 / 8
    @show min_ν = (1 / (1 / dx_min^2 + 1 / dy_min^2))^2 / timescale

    @inline function νhb_final(i, j, k, grid, clock, fields, p) 
        δ₁ = Dₛ(i, j, k, grid, fields.u, fields.v)    
        δ₂ = ℑxyᶜᶜᵃ(i, j, k, grid, Dₜ, fields.u, fields.v)    
        return max(Δᶜᶜᶜ(i, j, k, grid)^4 * p.C₄ * sqrt(δ₁^2 + δ₂^2), p.min_ν)
    end

    loc = (Center, Center, Center)

    return ScalarBiharmonicDiffusivity(formulation; 
                                       ν=νhb_final, discrete_form=true, loc, 
                                       parameters = (C₄ = C₄, min_ν = min_ν))
end

function leith_viscosity(formulation, grid; Cₗₜ = 1.0)

    dx_min    = min_Δx(grid.underlying_grid)
    dy_min    = min_Δy(grid.underlying_grid)
    timescale = 10days

    @show C₄    = (Cₗₜ / π)^3 / 8
    @show min_ν = (1 / (1 / dx_min^2 + 1 / dy_min^2))^2 / timescale

    @inline function νhb_final(i, j, k, grid, clock, fields, p) 
        ∂ζ = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)^2 + ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)^2
        ∂δ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)^2 + ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)^2

        return max(Δᶜᶜᶜ(i, j, k, grid)^5 * p.C₄ * sqrt(∂ζ + ∂δ), p.min_ν)
    end

    loc = (Center, Center, Center)

    return ScalarBiharmonicDiffusivity(formulation; 
                                       ν=νhb_final, discrete_form=true, loc, 
                                       parameters = (C₄ = C₄, min_ν = min_ν))
end