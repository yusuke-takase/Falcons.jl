@testset "orientation angle convention" begin
	# At (θ, φ) = (π/2, 0) the local East is +y and the local North is +z.
	@test orientation_angle(π / 2, 0.0, [0.0, 1.0, 0.0]) ≈ 0.0 atol = 1e-15
	@test orientation_angle(π / 2, 0.0, [0.0, 0.0, 1.0]) ≈ π / 2
	@test orientation_angle(π / 2, 0.0, [0.0, 0.0, -1.0]) ≈ -π / 2
	@test abs(orientation_angle(π / 2, 0.0, [0.0, -1.0, 0.0])) ≈ π
	# ψ is measured from East (e_φ) towards North (−e_θ).
	θ, φ = 1.1, 2.3
	e_φ = [-sin(φ), cos(φ), 0.0]
	e_θ = [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
	for ψ in range(-3.0, 3.0; length = 13)
		@test orientation_angle(θ, φ, cos(ψ) .* e_φ .- sin(ψ) .* e_θ) ≈ ψ atol = 1e-12
	end
end

@testset "compute_pointing_and_orientation" begin
	# litebird_sim's test: Ry(π/2) sends the beam to (π/2, 0) and x to −z (South).
	θ, φ, ψ = compute_pointing_and_orientation(quat_rotation_y(π / 2))
	@test [θ, φ, ψ] ≈ [π / 2, 0.0, -π / 2] atol = 1e-15
	θ, φ, ψ = compute_pointing_and_orientation(quat_rotation_x(π / 4) * quat_rotation_y(π / 2))
	@test ψ ≈ -π / 4
	for q in (quat_rotation_z(-0.3) * quat_rotation_y(1.2), quat_rotation_z(3.0) * quat_rotation_y(0.4))
		θ, φ, _ = compute_pointing_and_orientation(q)
		@test 0 <= φ < 2π
		@test [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)] ≈ rotate_vector(q, [0, 0, 1.0])
	end
end

@testset "spin to ecliptic" begin
	scan = SpinningScanningStrategy(spin_sun_angle_rad = deg2rad(45), precession_rate_hz = 1e-4, spin_rate_hz = 1e-3)
	# At t = 0 the spin axis lies in the x–z plane, α away from the anti-Sun direction x.
	@test rotate_vector(spin_to_ecliptic(scan, 0.0), [0, 0, 1.0]) ≈ [cos(deg2rad(45)), 0, sin(deg2rad(45))]
	# The spin axis is fixed by the spin rotation and the anti-Sun direction advances with Ω.
	t = 1000.0
	no_prec = SpinningScanningStrategy(spin_sun_angle_rad = deg2rad(45), precession_rate_hz = 0.0, spin_rate_hz = 1e-3)
	λ = Falcons.YEARLY_OMEGA_SPIN_HZ * t
	@test rotate_vector(spin_to_ecliptic(no_prec, t), [0, 0, 1.0]) ≈
		  rotate_vector(quat_rotation_z(λ), [cos(deg2rad(45)), 0, sin(deg2rad(45))])
	@test Falcons.YEARLY_OMEGA_SPIN_HZ == 2π / 31557600
end

@testset "satellite presets" begin
	scan, inst = satellite_preset(:LiteBIRD)
	@test scan.spin_sun_angle_rad ≈ deg2rad(45)
	@test scan.precession_rate_hz ≈ 1 / (192.348 * 60)
	@test scan.spin_rate_hz ≈ 0.05 / 60
	@test inst.spin_boresight_angle_rad ≈ deg2rad(50)
	@test :Planck in satellite_presets()
	@test_throws ArgumentError satellite_preset(:Unknown)
end

@testset "unit conversions" begin
	@test rpm2angfreq(60.0) ≈ 2π
	@test period2rpm(2.0) == 0.5
	@test period2rpm(30.0; unit = "sec") == 2.0
	@test rpm2period(0.5) == 2.0
	@test arcmin2rad(60.0) ≈ deg2rad(1.0)
end
