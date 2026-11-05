scan_obs, inst_obs = satellite_preset(:LiteBIRD)
dets_obs = [DetectorInfo(name = "a", sampling_rate_hz = 3.0), DetectorInfo(name = "b", sampling_rate_hz = 3.0,
	quat = detector_quat(0.05, 0.3))]

@testset "construction" begin
	obs = Observation(scan_obs, inst_obs, dets_obs; start_time_s = 10.0, duration_s = 4.5)
	@test obs.n_samples == 13  # floor(3 × 4.5)
	@test obs.sampling_rate_hz == 3.0
	@test n_detectors(obs) == 2
	@test get_times(obs) ≈ 10.0 .+ (0:12) ./ 3
	@test get_time_span(obs) ≈ 13 / 3
	@test obs.start_time_simulation_s == 10.0
	@test Observation(scan_obs, inst_obs, dets_obs[1]; start_time_s = 0.0, duration_s = 1.0).detectors == dets_obs[1:1]

	mixed = [DetectorInfo(name = "a", sampling_rate_hz = 3.0), DetectorInfo(name = "b", sampling_rate_hz = 4.0)]
	@test_throws ArgumentError Observation(scan_obs, inst_obs, mixed; start_time_s = 0.0, duration_s = 1.0)
	twins = [DetectorInfo(name = "a", sampling_rate_hz = 3.0), DetectorInfo(name = "a", sampling_rate_hz = 3.0)]
	@test_throws ArgumentError Observation(scan_obs, inst_obs, twins; start_time_s = 0.0, duration_s = 1.0)
	@test_throws ArgumentError Observation(scan_obs, inst_obs, dets_obs; start_time_s = 0.0, duration_s = 0.1)
	@test_throws ArgumentError Observation(scan_obs, inst_obs, dets_obs; start_time_s = 0.0, duration_s = 1.0, delta_time_s = 0.0)
end

@testset "split_observation" begin
	obs = Observation(scan_obs, inst_obs, dets_obs; start_time_s = 5.0, duration_s = 11 / 3,
		hwp = IdealHWP(ang_speed_radpsec = 4.8, start_angle_rad = 0.1))
	parts = split_observation(obs, 3)
	@test [p.n_samples for p in parts] == [4, 4, 3]
	@test [p.start_time_s for p in parts] ≈ [5.0, 5.0 + 4 / 3, 5.0 + 8 / 3]
	@test all(p.start_time_simulation_s == 5.0 for p in parts)
	@test reduce(vcat, get_times.(parts)) ≈ get_times(obs)
	whole, χ = get_pointings(obs)
	pieces = [get_pointings(p) for p in parts]
	@test reduce(vcat, last.(pieces)) ≈ χ atol = 1e-12
	@test reduce((a, b) -> cat(a, b; dims = 2), first.(pieces)) ≈ whole atol = 1e-12
	@test_throws ArgumentError split_observation(obs, 0)
	@test_throws ArgumentError split_observation(obs, 12)
end

@testset "IdealHWP" begin
	hwp = IdealHWP(ang_speed_radpsec = 2.0, start_angle_rad = 0.5)
	@test hwp_angle(hwp, 0.0) == 0.5
	@test hwp_angle(hwp, 1.0) ≈ rem2pi(2.5, RoundNearest)
	@test all(-π <= hwp_angle(hwp, t) <= π for t in 0:0.37:100)
	obs = Observation(scan_obs, inst_obs, dets_obs; start_time_s = 0.0, duration_s = 2.0)
	@test last(get_pointings(obs)) === nothing
end
