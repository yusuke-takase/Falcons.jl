"""
    Falcons

Fast pointing, hit-map and h-map (spin-moment) computation for spinning CMB
satellites, with the conventions, types and names of litebird_sim.

Main entry points:

- Scan description: [`SpinningScanningStrategy`](@ref), [`InstrumentInfo`](@ref),
  [`DetectorInfo`](@ref), [`IdealHWP`](@ref), [`Observation`](@ref),
  [`satellite_preset`](@ref), [`detector_quat`](@ref), [`split_observation`](@ref).
- Rotations: [`Quaternion`](@ref), [`RotQuaternion`](@ref) (constant or sampled
  in time, e.g. an external attitude), [`quat_rotation`](@ref).
- Pointings: [`get_pointings`](@ref), [`get_pointing_pixels`](@ref), [`get_hitmap`](@ref).
- h-maps: [`make_h_maps`](@ref), [`HMapsResult`](@ref), [`save_h_maps`](@ref),
  [`load_h_maps_from_file`](@ref).
- Instrument database: [`Imo`](@ref), [`query`](@ref), `DetectorInfo(imo, url)` and
  the other `(imo, url)` constructors, [`detector_list_from_parameters`](@ref).
- Parameter files: [`run_h_maps_pipeline`](@ref).
"""
module Falcons

using Base.Threads
using HDF5
using Healpix
import JSON
using LinearAlgebra
using Printf
using StaticArrays
import TOML

include("quaternions.jl")
include("coordinates.jl")
include("units.jl")
include("imo.jl")
include("scanning.jl")
include("detectors.jl")
include("hwp.jl")
include("observations.jl")
include("pointings.jl")
include("h_maps.jl")
include("pipelines.jl")

# quaternions.jl
export Quaternion, quat_rotation, quat_rotation_x, quat_rotation_y, quat_rotation_z, rotate_vector
export RotQuaternion, slerp, is_constant, time_range

# coordinates.jl
export CoordinateSystem, Ecliptic, Galactic
export rotate_coordinates_e2g!, rotate_vector_e2g!

# units.jl
export period2rpm, rpm2angfreq, rpm2period, arcmin2rad

# scanning.jl
export SpinningScanningStrategy, spin_to_ecliptic, orientation_angle, compute_pointing_and_orientation
export satellite_preset, satellite_presets

# imo.jl
export Imo, query, DataFile

# detectors.jl
export DetectorInfo, FreqChannelInfo, InstrumentInfo, bore2spin_quat, detector_quat
export get_boresight_detector, detector_list_from_parameters

# hwp.jl
export IdealHWP, hwp_angle

# observations.jl
export Observation, split_observation, get_times, get_time_span, n_detectors

# pointings.jl
export get_pointings, get_pointing_pixels, PointingPixels, get_hitmap

# h_maps.jl
export make_h_maps, HMapsResult, h_map, hitmap, nm_grid, save_h_maps, load_h_maps_from_file

# pipelines.jl
export run_h_maps_pipeline

end
