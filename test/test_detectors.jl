@testset "DetectorInfo" begin
	d = DetectorInfo()
	@test d.name == ""
	@test d.sampling_rate_hz == 0.0
	@test d.ellipticity == 1.0
	@test d.pol_efficiency == 1.0
	@test is_constant(d.quat)
	@test d.quat.quats[1] == one(Quaternion{Float64})

	q = detector_quat(0.1, 0.5; ψ = 0.2)
	@test DetectorInfo(name = "a", quat = q).quat.quats[1] ≈ q
	@test DetectorInfo(name = "a", quat = Falcons.components(q)).quat.quats[1] ≈ q
	series = RotQuaternion([q, q * quat_rotation_x(1e-3)]; start_time_s = 0.0, sampling_rate_hz = 1.0)
	@test !is_constant(DetectorInfo(name = "a", quat = series).quat)
	@test DetectorInfo(name = "a", sampling_rate_hz = 19).sampling_rate_hz === 19.0
	# Integer quaternion components.
	@test DetectorInfo(quat = [0, 0, 0, 1]).quat.quats[1] === one(Quaternion{Float64})
	@test DetectorInfo(quat = Quaternion(0, 0, 1, 0)).quat.quats[1] === Quaternion(0.0, 0.0, 1.0, 0.0)
	@test DetectorInfo(quat = Quaternion{Int}(0, 0, 0, 1)).quat.quats[1] === one(Quaternion{Float64})
end

@testset "detector_quat" begin
	θ, φ, ψ = 0.1, 0.5, 0.2
	q = detector_quat(θ, φ; ψ)
	@test q ≈ quat_rotation_z(φ) * quat_rotation_y(θ) * quat_rotation_z(ψ)
	# The beam points at (θ, φ) in the boresight frame.
	@test rotate_vector(q, [0, 0, 1.0]) ≈ [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]
end

@testset "InstrumentInfo" begin
	inst = InstrumentInfo(spin_boresight_angle_rad = 0.8, spin_rotangle_rad = 0.3, boresight_rotangle_rad = -0.4)
	@test bore2spin_quat(inst) ≈ quat_rotation_z(0.3) * quat_rotation_y(0.8) * quat_rotation_z(-0.4)
	@test bore2spin_quat(InstrumentInfo()) ≈ one(Quaternion{Float64})
	ch = FreqChannelInfo(bandcenter_ghz = 140.0)
	@test ch.fmin_hz == 1e-5
	@test ch.alpha == 1.0
end
