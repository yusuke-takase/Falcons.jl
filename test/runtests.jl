using Falcons
using HDF5
using Test

const REFDIR = joinpath(@__DIR__, "reference")

"Absolute difference between two angles, modulo 2π."
angdiff(a, b) = abs(rem2pi(a - b, RoundNearest))

@testset "Falcons.jl" begin
	@testset "quaternion" include("test_quaternion.jl")
	@testset "angles" include("test_angles.jl")
	@testset "invariants" include("test_invariants.jl")
	@testset "regression" include("test_regression.jl")
	@testset "io" include("test_io.jl")
	@testset "units" include("test_units.jl")
end
