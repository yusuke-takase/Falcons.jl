module Falcons

using ReferenceFrameRotations
using Healpix
using LinearAlgebra
using Base.Threads
using StaticArrays
using ProgressMeter

include("./function/scanning.jl")
<<<<<<< HEAD
=======
include("./function/mapmake.jl")
include("./function/make_scanmap.jl")
>>>>>>> 5f71a0a80c3b7b30cc36ace051d5bafac9ded6f1
include("./function/systematics.jl")

export ScanningStrategy, gen_ScanningStrategy, get_pointings, get_pointing_pixels, period2rpm, get_pointings_tuple, pixtod2hitmap, angtod2hitmap, convert_maps, rpm2angfreq
export get_psiDataBase, get_psi_time_DataBase, ScanningStrategy2MapInfo, set_input, orientation_func, get_hn_map, ScanningStrategy2map

end
