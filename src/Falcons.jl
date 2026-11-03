module Falcons

using Healpix
using LinearAlgebra
using Base.Threads
using StaticArrays
using JSON
using DataFrames
using Printf
using Statistics
using TOML
using HDF5

include("./function/quaternions.jl")
include("./function/scanning.jl")
include("./function/imo.jl")
include("./function/scanfields.jl")
include("./function/coordinates.jl")
include("./function/units.jl")
include("./function/pipelines.jl")

# quaternions.jl
export Quaternion, quat_rotation, quat_rotation_x, quat_rotation_y, quat_rotation_z, rotate_vector

# scanning.jl
export show_ss
export ScanningStrategy, pointings
export gen_ScanningStrategy, get_satellite
export get_pointings, get_pointing_pixels
export hitmap_from_angles, hitmap_from_pixels
export ang_to_quat, get_pointings_arbitrary_dets

# imo.jl
export Imo, gen_imo
export get_channel_list, get_channel_info, get_detectors
export imo_telescope!, imo_channel!, imo_name!, get_instrument_info
export get_pol_angle, imo2ecl_coordinates

# units.jl
export period2rpm, rpm2angfreq, rpm2period, arcmin2rad

# scanfields.jl
export ScanField, get_scanfield_step
export h_nm, get_hnm_quantify, get_scanfield

# coordinates.jl
export rotate_coordinates_e2g!, rotate_vector_e2g!

# pipelines.jl
export sim_det_scanfields, create_h5_file

end
