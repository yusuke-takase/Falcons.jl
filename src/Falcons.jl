module Falcons

using ReferenceFrameRotations
using Healpix
using LinearAlgebra
using Base.Threads
using StaticArrays

include("./function/func4scan.jl")
include("./function/scanning.jl")
include("./function/mapmake.jl")

export get_scan_tod, ScanStrategy
export Mapmaking, HitMap

end
