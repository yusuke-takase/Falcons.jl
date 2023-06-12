module Falcons

using ReferenceFrameRotations
using Healpix
using LinearAlgebra
using Base.Threads
using StaticArrays
using ProgressMeter
using JSON
using DataFrames
using Printf

include("./function/scanning.jl")
include("./function/systematics.jl")

export vec2ang_ver2, quaternion_rotator, imo2scan_coordinate, get_pointings_tuple, gen_ScanningStrategy_imo, ScanningStrategy_imo, imo_name!, imo_telescope!, imo_channel!
export ScanningStrategy, gen_ScanningStrategy, get_pointings, get_pointing_pixels, period2rpm, get_pointings_tuple, pixtod2hitmap, angtod2hitmap, convert_maps, rpm2angfreq, rpm2period, show_ss, ecliptic2galactic, pointings, normarize!
export get_psiDataBase, get_psi_time_DataBase, ScanningStrategy2MapInfo, set_input, orientation_func, get_hn_map, ScanningStrategy2map, binned_mapmake, signal, interp_signal

end
