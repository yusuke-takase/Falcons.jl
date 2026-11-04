@testset "durations" begin
	@test Falcons._duration_s(3600) == 3600.0
	@test Falcons._duration_s("1 year") == 365.25 * 86400
	@test Falcons._duration_s("2 days") == 2 * 86400
	@test Falcons._duration_s("30 min") == 1800
	@test Falcons._duration_s("1.5 hours") == 5400
	@test Falcons._duration_s("10 s") == 10
	@test_throws ArgumentError Falcons._duration_s("3 fortnights")
	@test_throws ArgumentError Falcons._start_time_s("2030-01-01")
end

@testset "run_h_maps_pipeline" begin
	imo = Imo(joinpath(REFDIR, "mock_imo"))
	rel = "/releases/vMOCK2/satellite"
	mktempdir() do dir
		params = joinpath(dir, "params.toml")
		write(params, """
		[simulation]
		base_path = "$(escape_string(dir))"
		start_time = 100.0
		duration_s = "1 min"

		[imo]
		flatfile_location = "$(escape_string(joinpath(REFDIR, "mock_imo")))"

		[scanning_strategy]
		imo_url = "$rel/scanning_parameters"

		[instrument]
		imo_url = "$rel/MOCKFT/instrument_info"

		[hwp]
		ang_speed_radpsec = 4.8

		[[detectors]]
		channel_info_obj = "$rel/MOCKFT/M-100/channel_info"
		num_of_detectors_from_channel = 2

		[[detectors]]
		name = "extra"
		sampling_rate_hz = 19.0

		[h_maps]
		nside = 8
		n_m_couples = [[0, 0], [2, 0], [2, -4]]
		output_coordinate_system = "Ecliptic"
		output_directory = "maps"
		""")
		res = run_h_maps_pipeline(params)
		@test res.detectors == ["M100_00_A", "M100_00_B", "extra"]
		@test res.coordinate_system == Ecliptic
		@test sort(readdir(joinpath(dir, "maps"))) == ["h_maps_det_M100_00_A.h5", "h_maps_det_M100_00_B.h5", "h_maps_det_extra.h5"]

		dets = [DetectorInfo(imo, "$rel/MOCKFT/M-100/M100_00_A/detector_info"), DetectorInfo(imo, "$rel/MOCKFT/M-100/M100_00_B/detector_info"),
			DetectorInfo(name = "extra", sampling_rate_hz = 19.0)]
		obs = Observation(SpinningScanningStrategy(imo, "$rel/scanning_parameters"), InstrumentInfo(imo, "$rel/MOCKFT/instrument_info"),
			dets; start_time_s = 100.0, duration_s = 60.0, hwp = IdealHWP(ang_speed_radpsec = 4.8))
		direct = make_h_maps(obs; nside = 8, n_m_couples = [(0, 0), (2, 0), (2, -4)], output_coordinate_system = Ecliptic)
		@test res.h == direct.h
		back = load_h_maps_from_file(joinpath(dir, "maps", "h_maps_det_extra.h5"))
		@test h_map(back, "extra", 2, -4) == h_map(direct, "extra", 2, -4)

		write(params, """
		[simulation]
		duration_s = 60
		[scanning_strategy]
		imo_url = "$rel/scanning_parameters"
		[[detectors]]
		name = "d"
		sampling_rate_hz = 1.0
		[h_maps]
		nside = 8
		""")
		@test_throws ArgumentError run_h_maps_pipeline(params)  # imo_url without [imo]
	end
end
