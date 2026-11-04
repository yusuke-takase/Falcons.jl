using Falcons
using HDF5
using Healpix
using LinearAlgebra
using StaticArrays
using Test

const REFDIR = joinpath(@__DIR__, "reference")

"Absolute difference between two angles, modulo 2π."
angdiff(a, b) = abs(rem2pi(a - b, RoundNearest))

include("lbs_reference.jl")

@testset "Falcons.jl" begin
	@testset "quaternions" include("test_quaternions.jl")
	@testset "units" include("test_units.jl")
	@testset "scanning" include("test_scanning.jl")
	@testset "detectors" include("test_detectors.jl")
	@testset "coordinates" include("test_coordinates.jl")
	@testset "healpix" include("test_healpix.jl")
	@testset "observations" include("test_observations.jl")
	@testset "pointings" include("test_pointings.jl")
	@testset "h_maps" include("test_h_maps.jl")
	@testset "imo" include("test_imo.jl")
	@testset "pipelines" include("test_pipelines.jl")
end
