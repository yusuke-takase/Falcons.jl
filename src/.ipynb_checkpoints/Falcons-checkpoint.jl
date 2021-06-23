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

export ScanningStrategy, gen_ScanningStrategy, get_pointings, get_pointing_pixels, period2rpm
export TwoTelescopes_ScanningStrategy2map, ThreeTelescopes_ScanningStrategy2map
export ScanningStrategy2map, Mapmaking, pixtod2hitmap, angtod2hitmap, Genmap

end
