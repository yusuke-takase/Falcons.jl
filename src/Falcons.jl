module Falcons

using ReferenceFrameRotations
using Healpix
using LinearAlgebra
using Base.Threads
using StaticArrays
using ProgressMeter
using Printf

include("./function/scanning.jl")
include("./function/systematics.jl")

export ScanningStrategy, gen_ScanningStrategy, get_pointings, get_pointing_pixels, period2rpm, get_pointings_tuple, pixtod2hitmap, angtod2hitmap, convert_maps, rpm2angfreq, rpm2period, show_ss, ecliptic2galactic, pointings, normarize!
export get_psiDataBase, get_psi_time_DataBase, ScanningStrategy2MapInfo, set_input, orientation_func, get_hn_map, ScanningStrategy2map, binned_mapmake, signal, interp_signal

end
