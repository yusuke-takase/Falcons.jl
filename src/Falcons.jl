module Falcons

using ReferenceFrameRotations
using Healpix
using LinearAlgebra
using Base.Threads
using StaticArrays
using ProgressMeter

include("./function/func4scan.jl")
include("./function/scanning.jl")
include("./function/systematics.jl")
include("./function/mapmake.jl")
include("./function/make_scanmap.jl")

#include("./function/focalplane.jl")

export ScanningStrategy, gen_ScanningStrategy, get_pointings, get_pointing_pixels, period2rpm, get_pointings_xyz_tuple, get_pointings_tuple
export TwoTelescopes_ScanningStrategy2map, TwoTelescopes_ScanningStrategy2MapInfo, ThreeTelescopes_ScanningStrategy2map
export ScanningStrategy2map, Mapmaking, pixtod2hitmap, angtod2hitmap, array2map, xyztod2hitmap, convert_maps, get_hn_map
export get_psiDataBase, get_psi_time_DataBase, ScanningStrategy2MapInfo, rpm2angfreq
export set_input, orientation_func, get_hmap
#export w_μ, H_μ, hit_matrix

end
