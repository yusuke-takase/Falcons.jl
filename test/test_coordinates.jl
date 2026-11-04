@testset "Ecliptic to Galactic" begin
	M = Falcons.ECL_TO_GAL_ROT_MATRIX
	@test M * M' ≈ I atol = 1e-14
	# The quaternion and the matrix describe the same rotation.
	for v in ([1.0, 0, 0], [0, 1.0, 0], [0.3, -0.4, 0.866])
		@test rotate_vector(Falcons.ECL_TO_GAL_QUAT, v) ≈ M * v atol = 1e-14
	end
	# The north ecliptic pole sits at (l, b) ≈ (96.3843°, 29.8113°).
	v = rotate_vector_e2g!([0.0, 0.0, 1.0])
	@test v ≈ Falcons.NORTH_POLE_VEC
	@test rad2deg(atan(v[2], v[1])) ≈ 96.3843 atol = 1e-3
	@test rad2deg(asin(v[3])) ≈ 29.8113 atol = 1e-3
	# Same values as litebird_sim.
	h5open(LBS_REF) do f
		@test Falcons.components(Falcons.ECL_TO_GAL_QUAT) == read(f["config/ecl_to_gal_quat"])
		@test M == permutedims(read(f["config/ecl_to_gal_rot_matrix"]))
	end
end

@testset "rotate_coordinates_e2g! agrees with Galactic pointings" begin
	scan, inst = satellite_preset(:LiteBIRD)
	det = DetectorInfo(name = "d", sampling_rate_hz = 1.0, quat = detector_quat(0.05, 1.0; ψ = 0.4))
	obs = Observation(scan, inst, det; start_time_s = 0.0, duration_s = 600.0)
	ecl, _ = get_pointings(obs; output_coordinate_system = Ecliptic)
	gal, _ = get_pointings(obs; output_coordinate_system = Galactic)
	θ, φ, ψ = ecl[1, :, 1], ecl[2, :, 1], ecl[3, :, 1]
	rotate_coordinates_e2g!(θ, φ, ψ)
	@test maximum(abs.(θ .- gal[1, :, 1])) < 1e-10
	@test maximum(angdiff.(φ, gal[2, :, 1])) < 1e-10
	@test maximum(angdiff.(ψ, gal[3, :, 1])) < 1e-10
end
