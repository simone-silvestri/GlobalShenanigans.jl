module GlobalShenanigans

using Reexport
@reexport using Oceananigans
using JLD2
using Oceananigans.Fields: interpolate
using KernelAbstractions: @kernel, @index, MultiEvent
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device, arch_array, device_event
using Oceananigans.Units
using Oceananigans.Fields: interpolate
using Statistics: dot
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Oceananigans.Grids: architecture

const ABOVE_SEA_LEVEL = 100

include("src/utils.jl")
include("src/interpolate_bathymetry.jl")
include("src/interpolate_initial_condition.jl")
# include("src/interpolate_fluxes.jl")

end # module
