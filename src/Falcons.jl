module Falcons

using ReferenceFrameRotations
using Healpix
using LinearAlgebra
using Base.Threads
using StaticArrays
using ProgressMeter
using JSON
using DataFrames
using SkyCoords
using Printf
using Statistics
using TOML
using HDF5
using Base.Filesystem
using NaNStatistics

include("./function/scanning.jl")
include("./function/imo.jl")
include("./function/scanfields.jl")
include("./function/coordinates.jl")
include("./function/units.jl")
include("./function/pipelines.jl")

# scanning.jl
export rotate_quat, show_ss
export ScanningStrategy, pointings
export gen_ScanningStrategy, get_satellite
export get_pointings, get_pointing_pixels
export angtod2hitmap, pixtod2hitmap, normarize!
export ang_to_quat, get_pointings_arbitrary_dets

# imo.jl
export Imo, gen_imo
export get_channel_list, get_channel_info, get_detectors
export imo_telescope!, imo_channel!, imo_name!, get_instrument_info
export get_pol_angle, imo2ecl_coordinates

# units.jl
export period2rpm, rpm2angfreq, rpm2period, arcmin2rad
export convert_maps

# scanfields.jl
export scanfield, get_scanfield_step
export h_nm, get_hnm_quantify, get_scanfield

# coordinates.jl
export _ang2galvec_one_sample, rotate_coordinates_e2g!, rotate_vec_ecl_to_gal_coordinates!, _rotate_coordinates_and_pol_e2g_for_one_sample

# pipelines.jl
export sim_det_scanfields, create_h5_file

end
