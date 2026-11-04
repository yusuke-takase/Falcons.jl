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
	@test abs(normalize(q)) ≈ 1
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
	axis = SVector(0.3, -0.5, 0.8) / norm(SVector(0.3, -0.5, 0.8))
	v = SVector(0.7, 0.1, -0.2)
	θ = 1.234
	rodrigues = v * cos(θ) + cross(axis, v) * sin(θ) + axis * dot(axis, v) * (1 - cos(θ))
	q = quat_rotation(θ, axis)
	@test rotate_vector(q, v) ≈ rodrigues
	@test Falcons.rotate_x_vector(q) ≈ rotate_vector(q, ex)
	@test Falcons.rotate_z_vector(q) ≈ rotate_vector(q, ez)
	# (q1 * q2) rotates by q2 first, then q1.
	q1, q2 = quat_rotation_z(0.4), quat_rotation_x(-1.1)
	@test rotate_vector(q1 * q2, v) ≈ rotate_vector(q1, rotate_vector(q2, v))
end

@testset "slerp" begin
	a, b = quat_rotation_z(0.2), quat_rotation_z(1.4)
	@test slerp(a, b, 0.0) ≈ a
	@test slerp(a, b, 1.0) ≈ b
	@test slerp(a, b, 0.25) ≈ quat_rotation_z(0.5)
	# Shortest path: −b is the same rotation.
	@test slerp(a, -b, 0.25) ≈ quat_rotation_z(0.5)
	@test slerp(a, a, 0.3) ≈ a
end

quats_approx(a, b) = length(a) == length(b) && all(isapprox.(a, b))

@testset "RotQuaternion" begin
	c = RotQuaternion(quat_rotation_y(0.3))
	@test is_constant(c)
	@test slerp(c, 1e9) == c.quats[1]
	@test RotQuaternion([0.0, 0.0, 0.0, 2.0]).quats[1] == one(Quaternion{Float64})  # normalised
	@test_throws ArgumentError RotQuaternion(Quaternion{Float64}[]; start_time_s = 0.0, sampling_rate_hz = 1.0)

	quats = [quat_rotation_z(0.1 * k) for k in 0:10]
	r = RotQuaternion(quats; start_time_s = 100.0, sampling_rate_hz = 2.0)
	@test !is_constant(r)
	@test time_range(r) == (100.0, 105.0)
	@test slerp(r, 100.0) ≈ quats[1]
	@test slerp(r, 105.0) ≈ quats[end]
	@test slerp(r, 101.25) ≈ quat_rotation_z(0.25)
	@test_throws ArgumentError slerp(r, 99.0)
	@test_throws ArgumentError slerp(r, 105.5)
	@test_throws ArgumentError RotQuaternion(quats[1:1]; start_time_s = 0.0, sampling_rate_hz = 1.0)
	@test_throws ArgumentError RotQuaternion(quats; start_time_s = 0.0, sampling_rate_hz = 0.0)

	m = reduce(hcat, Falcons.components.(quats))
	@test quats_approx(RotQuaternion(m; start_time_s = 100.0, sampling_rate_hz = 2.0).quats, r.quats)
	@test_throws ArgumentError RotQuaternion(permutedims(m); start_time_s = 0.0, sampling_rate_hz = 1.0)

	g = quat_rotation_x(0.7)
	@test quats_approx((g * r).quats, [g * q for q in quats])
	@test quats_approx((r * g).quats, [q * g for q in quats])
	@test quats_approx((c * r).quats, [c.quats[1] * q for q in quats])
	@test quats_approx((r * r).quats, [q * q for q in quats])
	other = RotQuaternion(quats; start_time_s = 101.0, sampling_rate_hz = 2.0)
	@test_throws ArgumentError r * other
end
