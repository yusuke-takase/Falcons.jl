using LinearAlgebra: cross, dot, norm
using StaticArrays: SVector

@testset "Hamilton product" begin
	i = Quaternion(1.0, 0.0, 0.0, 0.0)
	j = Quaternion(0.0, 1.0, 0.0, 0.0)
	k = Quaternion(0.0, 0.0, 1.0, 0.0)
	minus_one = Quaternion(0.0, 0.0, 0.0, -1.0)
	@test i * j == k
	@test j * k == i
	@test k * i == j
	@test i * i == minus_one
	@test i * j * k == minus_one
	@test j * i == Quaternion(0.0, 0.0, -1.0, 0.0)
end

@testset "algebra" begin
	q = Quaternion(0.1, -0.2, 0.3, 0.9)
	r = Quaternion(-0.4, 0.5, 0.2, -0.1)
	@test Quaternion([0.1, -0.2, 0.3, 0.9]) == q
	@test_throws ArgumentError Quaternion([1.0, 2.0, 3.0])
	@test Falcons.components(q) == SVector(0.1, -0.2, 0.3, 0.9)
	@test q * one(q) == q
	@test one(q) * q == q
	@test q * inv(q) ≈ one(q)
	@test conj(q * r) ≈ conj(r) * conj(q)
	@test abs(q * r) ≈ abs(q) * abs(r)
	@test isbitstype(Quaternion{Float64})
end

@testset "active right-handed rotations" begin
	ex, ey, ez = SVector(1.0, 0, 0), SVector(0, 1.0, 0), SVector(0, 0, 1.0)
	@test rotate_vector(quat_rotation(π / 2, ez), ex) ≈ ey
	@test rotate_vector(quat_rotation(π / 2, ex), ey) ≈ ez
	@test rotate_vector(quat_rotation(π / 2, ey), ez) ≈ ex
	for θ in (-2.1, 0.3, 1.7)
		@test quat_rotation_x(θ) ≈ quat_rotation(θ, ex)
		@test quat_rotation_y(θ) ≈ quat_rotation(θ, ey)
		@test quat_rotation_z(θ) ≈ quat_rotation(θ, ez)
	end
	# Rodrigues' formula for a generic axis and angle.
	axis = SVector(0.3, -0.5, 0.8) / norm(SVector(0.3, -0.5, 0.8))
	v = SVector(0.7, 0.1, -0.2)
	θ = 1.234
	rodrigues = v * cos(θ) + cross(axis, v) * sin(θ) + axis * dot(axis, v) * (1 - cos(θ))
	@test rotate_vector(quat_rotation(θ, axis), v) ≈ rodrigues
	# Composition: (q1 * q2) rotates by q2 first, then q1.
	q1, q2 = quat_rotation_z(0.4), quat_rotation_x(-1.1)
	@test rotate_vector(q1 * q2, v) ≈ rotate_vector(q1, rotate_vector(q2, v))
end
