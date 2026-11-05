const LBS_ATOL_ANGLE = 1e-10 # rad
const LBS_ATOL_HWP = 1e-11   # rad; litebird_sim accumulates the HWP angle with a recurrence

function check_against_lbs(f, case, obs; hwp = true)
	g = f[case]
	idx = read(g["sample_index"])
	@test obs.n_samples == read(attributes(g)["n_samples"])
	@test obs.start_time_s == read(attributes(g)["start_time_s"])
	pt_ecl, χ = get_pointings(obs; output_coordinate_system = Ecliptic)
	pt_gal, _ = get_pointings(obs; output_coordinate_system = Galactic)
	for (pt, key) in ((pt_ecl, "pointings_ecl"), (pt_gal, "pointings_gal"))
		dθ, dφ, dψ = lbs_pointing_errors(pt, read(g[key]), idx)
		@test dθ < LBS_ATOL_ANGLE
		@test dφ < LBS_ATOL_ANGLE
		@test dψ < LBS_ATOL_ANGLE
	end
	if hwp
		@test maximum(angdiff.(χ[idx .+ 1], read(g["hwp_angle"]))) < LBS_ATOL_HWP
	end
end

h5open(LBS_REF) do f
	fs = lbs_cfg(f, "sampling_rate_hz")
	T = lbs_cfg(f, "duration_s")
	t0 = lbs_cfg(f, "start_time_b")

	@testset "A: boresight, start at 0" begin
		obs = Observation(lbs_scan(f), lbs_instrument_a(f), DetectorInfo(name = "boresight", sampling_rate_hz = fs);
			start_time_s = 0.0, duration_s = T)
		check_against_lbs(f, "A", obs; hwp = false)
		pt, _ = get_pointings(obs)
		@test pt[:, 1, 1] ≈ [deg2rad(95.0), 0.0, -π / 2]
	end

	@testset "B: rotated instrument, offset detectors, HWP" begin
		check_against_lbs(f, "B", lbs_observation_b(f))
	end

	@testset "C: second of two observations" begin
		parts = split_observation(lbs_observation_b(f), 2)
		@test parts[2].start_time_simulation_s == read(attributes(f["C"])["start_time_simulation_s"])
		check_against_lbs(f, "C", parts[2])
	end

	@testset "D: interpolated pointings (Δ = 2 s)" begin
		check_against_lbs(f, "D", lbs_observation_b(f; delta_time_s = 2.0))
	end

	g = f["G_exact"]
	det2_quats = RotQuaternion(read(g["det2_quats"]); start_time_s = read(attributes(g)["det2_start_time_s"]),
		sampling_rate_hz = read(attributes(g)["det2_sampling_rate_hz"]))
	dets_g = lbs_detectors_b(f)
	dets_g[3] = DetectorInfo(name = "det2", sampling_rate_hz = fs, pol_angle_rad = dets_g[3].pol_angle_rad, quat = det2_quats)
	obs_g(; kw...) = Observation(lbs_scan(f), lbs_instrument_b(f), dets_g; start_time_s = t0, duration_s = T,
		hwp = lbs_hwp_b(f), kw...)

	@testset "G: time-dependent detector quaternion" begin
		check_against_lbs(f, "G_exact", obs_g(); hwp = false)
		check_against_lbs(f, "G_interp", obs_g(delta_time_s = 2.0); hwp = false)
	end

	@testset "H: attitude given as a time series" begin
		h = f["H"]
		attitude = RotQuaternion(read(h["attitude_quats"]); start_time_s = read(attributes(h)["attitude_start_time_s"]),
			sampling_rate_hz = read(attributes(h)["attitude_sampling_rate_hz"]))
		obs = Observation(attitude, lbs_instrument_b(f), lbs_detectors_b(f); start_time_s = t0, duration_s = T)
		check_against_lbs(f, "H", obs; hwp = false)
	end

	@testset "I: spin-frame perturbation" begin
		i = f["I"]
		perturbation = RotQuaternion(read(i["perturbation_quats"]); start_time_s = read(attributes(i)["perturbation_start_time_s"]),
			sampling_rate_hz = read(attributes(i)["perturbation_sampling_rate_hz"]))
		check_against_lbs(f, "I", lbs_observation_b(f; spin_perturbation = perturbation); hwp = false)
	end

	@testset "F: hit map" begin
		@test get_hitmap(lbs_observation_b(f); nside = 32) == read(f["F/hitmap"])
	end

	@testset "time coverage is checked" begin
		short = RotQuaternion(read(f["H/attitude_quats"])[:, 1:100]; start_time_s = t0, sampling_rate_hz = 1.0)
		obs = Observation(short, lbs_instrument_b(f), lbs_detectors_b(f); start_time_s = t0, duration_s = T)
		@test_throws ArgumentError get_pointings(obs)
		late = RotQuaternion(read(f["G_exact/det2_quats"]); start_time_s = t0 + 10.0, sampling_rate_hz = 1.0)
		dets = [DetectorInfo(name = "x", sampling_rate_hz = fs, quat = late)]
		@test_throws ArgumentError get_pointings(Observation(lbs_scan(f), lbs_instrument_b(f), dets; start_time_s = t0, duration_s = T))
	end
end

@testset "detector selection" begin
	scan, inst = satellite_preset(:LiteBIRD)
	dets = [DetectorInfo(name = "a", sampling_rate_hz = 2.0), DetectorInfo(name = "b", sampling_rate_hz = 2.0, quat = detector_quat(0.1, 0.2))]
	obs = Observation(scan, inst, dets; start_time_s = 0.0, duration_s = 50.0)
	all_pt, _ = get_pointings(obs)
	@test size(all_pt) == (3, 100, 2)
	@test get_pointings(obs; detector_idx = 2)[1] == all_pt[:, :, 2]
	@test get_pointings(obs; detector_idx = [2, 1])[1] == all_pt[:, :, [2, 1]]
	@test_throws ArgumentError get_pointings(obs; detector_idx = 3)
end

@testset "pointing pixels and hit maps" begin
	scan, inst = satellite_preset(:LiteBIRD)
	dets = [DetectorInfo(name = "a", sampling_rate_hz = 5.0), DetectorInfo(name = "b", sampling_rate_hz = 5.0, quat = detector_quat(0.1, 0.2))]
	obs = Observation(scan, inst, dets; start_time_s = 3.0, duration_s = 3600.0, hwp = IdealHWP(ang_speed_radpsec = 1.0))
	pt, χ = get_pointings(obs; output_coordinate_system = Galactic)
	for nest in (false, true)
		px = get_pointing_pixels(obs; nside = 64, nest)
		resol = Resolution(64)
		expected = [nest ? ang2pixNest(resol, pt[1, k, d], pt[2, k, d]) : ang2pixRing(resol, pt[1, k, d], pt[2, k, d])
					for k in 1:obs.n_samples, d in 1:2]
		@test px.pixels == expected
		@test px.psi == permutedims(pt[3, :, :], (1, 2))
		@test px.hwp_angle == χ
		@test (px.nside, px.nest, px.coordinate_system, px.detectors) == (64, nest, Galactic, ["a", "b"])
		hits = get_hitmap(obs; nside = 64, nest)
		@test sum(hits) == 2 * obs.n_samples
		@test hits == [count(==(p), px.pixels) for p in 1:nside2npix(64)]
	end
	@test get_hitmap(split_observation(obs, 3); nside = 64) == get_hitmap(obs; nside = 64)
end

@testset "interpolation anchors are computed per block" begin
	scan, inst = satellite_preset(:LiteBIRD)
	q0 = detector_quat(0.1, 0.2)
	wobble = RotQuaternion([quat_rotation_x(1e-4 * sin(0.3 * k)) * q0 for k in 0:40]; start_time_s = 99.0, sampling_rate_hz = 0.25)
	dets = [DetectorInfo(name = "a", sampling_rate_hz = 3.0), DetectorInfo(name = "b", sampling_rate_hz = 3.0, quat = q0),
		DetectorInfo(name = "c", sampling_rate_hz = 3.0, quat = wobble)]
	obs = Observation(scan, inst, dets; start_time_s = 100.0, duration_s = 150.0, hwp = IdealHWP(ang_speed_radpsec = 1.0),
		start_time_simulation_s = 7.0, delta_time_s = 2.0)
	ref_pt, _ = get_pointings(obs)
	ref_px = get_pointing_pixels(obs; nside = 32)
	ref_hits = get_hitmap(obs; nside = 32)
	ref_h = make_h_maps(obs; nside = 16, n_m_couples = [(0, 0), (2, 0), (2, 4)])
	# The plan holds no anchors, whatever the duration of the observation.
	long = Observation(scan, inst, dets[1:2]; start_time_s = 100.0, duration_s = 30 * 86400.0, delta_time_s = 2.0)
	short = Observation(scan, inst, dets[1:2]; start_time_s = 100.0, duration_s = 60.0, delta_time_s = 2.0)
	@test Base.summarysize(Falcons._pointing_plan(long, Galactic, :all)) ==
		  Base.summarysize(Falcons._pointing_plan(short, Galactic, :all))
	default = Falcons._CHUNK_SAMPLES[]
	try
		for len in (1, 2, 7, 64)
			Falcons._CHUNK_SAMPLES[] = len
			@test get_pointings(obs)[1] == ref_pt
			px = get_pointing_pixels(obs; nside = 32)
			@test (px.pixels, px.psi) == (ref_px.pixels, ref_px.psi)
			@test get_hitmap(obs; nside = 32) == ref_hits
			h = make_h_maps(obs; nside = 16, n_m_couples = [(0, 0), (2, 0), (2, 4)])
			@test h.hits == ref_h.hits
			@test h.h == ref_h.h
			blk = Falcons._anchor_block(Falcons._pointing_plan(obs, Galactic, :all), 1, len)
			@test size(blk.quats, 1) <= cld(len, 6) + 2  # 6 samples per anchor interval
			# Samples sparser than the anchors (one sample every 4 anchors).
			sparse = Observation(scan, inst, DetectorInfo(name = "s", sampling_rate_hz = 0.125); start_time_s = 100.0,
				duration_s = 4000.0, delta_time_s = 2.0)
			plan = Falcons._pointing_plan(sparse, Galactic, :all)
			@test size(Falcons._anchor_block(plan, 1, Falcons._block_length(plan)).quats, 1) <= len + 4
		end
	finally
		Falcons._CHUNK_SAMPLES[] = default
	end
end
