module Falcons

using ReferenceFrameRotations
using Healpix
using LinearAlgebra
using Base.Threads
using StaticArrays
using ProgressMeter

include("./function/func4scan.jl")
include("./function/scanning.jl")
include("./function/mapmake.jl")
include("./function/make_scanmap.jl")
#include("./function/focalplane.jl")

export ScanningStrategy, gen_ScanningStrategy, get_pointings, get_pointing_pixels, period2rpm
export TwoTelescopes_ScanningStrategy2map, TwoTelescopes_ScanningStrategy2MapInfo, ThreeTelescopes_ScanningStrategy2map
export ScanningStrategy2map, Mapmaking, pixtod2hitmap, angtod2hitmap, Genmap
export get_psiDataBase, ScanningStrategy2MapInfo, get_pointings_tuple
#export lft_focalplane_configration, mft_focalplane_configration, hft_focalplane_configration, get_FP_each_freq
#export pickup_wafer, pickup_freq, ang2xy, longitude2pix, ang2xy_projected_hft, ang2xy_projected_mft

end
