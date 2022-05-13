using Oceananigans
using JLD2
using Plots
using Oceananigans.Fields: interpolate
using KernelAbstractions: @kernel, @index, MultiEvent
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device, arch_array, device_event

const ABOVE_SEA_LEVEL =  100
const degree   = 1
const latitude = 75

include("utils.jl")
include("interpolate_bathymetry.jl")

include("interpolate_fluxes.jl")
include("interpolate_initial_condition.jl")

## Start of the interpolation

arch = GPU()
bathymetry_original = interpolate_bathymetry_from_file("bathymetry-ice-21600x10800.jld2", 15)

## Start the modifications

bathy = Vector(undef, 10)

# mask individual elements

bathy[1]  = mask_isolated_elements(bathymetry_original, arch)
bathy[2]  = mask_coupled_elements(bathy[1], arch)
bathy[3]  = fill_up_all_lakes(bathy[2], arch, 40, 40)
bathy[4]  = mask_isolated_elements(bathy[3], arch)
bathy[5]  = mask_coupled_elements(bathy[4], arch)
bathy[6]  = fill_up_all_lakes(bathy[5], arch, 50, 50)
bathy[7]  = mask_isolated_elements(bathy[6], arch)
bathy[8]  = mask_coupled_elements(bathy[7], arch)
bathy[9]  = fill_up_all_lakes(bathy[8], arch, 50, 50)
bathy[10] = mask_isolated_elements(bathy[9], arch)
