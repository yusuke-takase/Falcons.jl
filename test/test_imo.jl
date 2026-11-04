const MOCK_IMO = joinpath(REFDIR, "mock_imo")
const MOCK_REL = "/releases/vMOCK2/satellite"

imo = Imo(MOCK_IMO)

@testset "loading" begin
	@test Imo(joinpath(MOCK_IMO, "schema.json")).location == imo.location
	@test sort(collect(keys(imo.releases))) == ["vMOCK1", "vMOCK2"]
	mktempdir() do dir
		@test_throws ArgumentError Imo(dir)
		touch(joinpath(dir, "schema.json.gz"))
		err = try
			Imo(dir)
		catch e
			e
		end
		@test err isa ArgumentError && occursin("decompress", err.msg)
	end
	@test_throws ArgumentError Imo(joinpath(MOCK_IMO, "missing.json"))
end

@testset "query" begin
	url = "$MOCK_REL/MOCKFT/M-100/M100_00_B/detector_info"
	df = query(imo, url)
	@test df isa DataFile
	@test df.metadata["name"] == "M100_00_B"
	@test query(imo, replace(url, "/releases" => ""; count = 1)) === df  # without the /releases prefix
	@test query(imo, url * "/").uuid == df.uuid                       # trailing slash
	@test query(imo, df.uuid) === df
	@test query(imo, uppercase(df.uuid)) === df
	@test query(imo, "/data_files/$(df.uuid)") === df
	@test query(imo, "/quantities/$(df.quantity)") isa Falcons.Quantity
	@test query(imo, "/entities/$(imo.path_to_entity["/satellite/MOCKFT"])").name == "MOCKFT"
	@test sort(df.release_tags) == ["vMOCK1", "vMOCK2"]
	# Release selection: vMOCK1 holds an older version of this detector.
	old = query(imo, "/releases/vMOCK1/satellite/MOCKFT/M-100/M100_00_A/detector_info")
	new = query(imo, "$MOCK_REL/MOCKFT/M-100/M100_00_A/detector_info")
	@test old.uuid != new.uuid
	@test old.release_tags == ["vMOCK1"]
	@test_throws ArgumentError query(imo, "/releases/vNONE/satellite/scanning_parameters")
	@test_throws ArgumentError query(imo, "$MOCK_REL/NOPE/instrument_info")
	@test_throws ArgumentError query(imo, "$MOCK_REL/MOCKFT/no_such_quantity")
	@test_throws ArgumentError query(imo, "/releases/vMOCK2")
	@test_throws ArgumentError query(imo, "00000000-0000-0000-0000-000000000000")
end

h5open(LBS_REF) do f
	g = f["imo"]
	ga = attributes(g)

	@testset "InstrumentInfo and SpinningScanningStrategy against litebird_sim" begin
		inst = InstrumentInfo(imo, "$MOCK_REL/MOCKFT/instrument_info")
		for key in ("boresight_rotangle_rad", "spin_boresight_angle_rad", "spin_rotangle_rad", "hwp_rpm", "wafer_space_cm")
			@test getfield(inst, Symbol(key)) ≈ read(ga["instrument_$key"]) rtol = 1e-15
		end
		@test Falcons.components(bore2spin_quat(inst)) ≈ read(g["instrument_bore2spin_quat"]) atol = 1e-15
		@test inst.name == "MOCKFT"
		@test inst.wafer_names == ["W0", "W1"]
		@test length(inst.channel_objs) == 2
		scan = SpinningScanningStrategy(imo, "$MOCK_REL/scanning_parameters")
		for key in ("spin_sun_angle_rad", "precession_rate_hz", "spin_rate_hz")
			@test getfield(scan, Symbol(key)) ≈ read(ga["scan_$key"]) rtol = 1e-15
		end
	end

	@testset "FreqChannelInfo and DetectorInfo against litebird_sim" begin
		for channel in ("M-100", "M-200")
			c = g["channel/$channel"]
			ca = attributes(c)
			ch = FreqChannelInfo(imo, "$MOCK_REL/MOCKFT/$channel/channel_info")
			@test ch.channel == channel
			for key in ("number_of_detectors", "net_detector_ukrts", "net_channel_ukrts", "bandcenter_ghz", "fmin_hz", "alpha")
				@test getfield(ch, Symbol(key)) ≈ read(ca[key]) rtol = 1e-15
			end
			@test ch.detector_objs == read(ca["detector_objs"])
			for name in keys(c)
				d = c[name]
				da = attributes(d)
				det = DetectorInfo(imo, "$MOCK_REL/MOCKFT/$channel/$name/detector_info")
				@test det.name == name
				for key in ("sampling_rate_hz", "fwhm_arcmin", "bandcenter_ghz", "net_ukrts", "fknee_mhz", "pol_angle_rad", "psi_rad", "squid")
					@test getfield(det, Symbol(key)) == read(da[key])
				end
				for key in ("wafer", "pixel", "pol", "orient", "channel")
					@test string(getfield(det, Symbol(key))) == read(da[key])
				end
				quats = read(d["quats"])  # (4, n)
				@test reduce(hcat, Falcons.components.(det.quat.quats)) ≈ quats atol = 1e-15
				if haskey(da, "quat_start_time_s")
					@test det.quat.start_time_s == read(da["quat_start_time_s"])
					@test det.quat.sampling_rate_hz == read(da["quat_sampling_rate_hz"])
				else
					@test is_constant(det.quat)
				end
			end
		end
	end

	@testset "J: observation described by the IMO" begin
		j = f["J"]
		dets = detector_list_from_parameters(imo, [Dict("channel_info_obj" => "$MOCK_REL/MOCKFT/M-100/channel_info")])
		@test [d.name for d in dets] == read(attributes(j)["detector_names"])
		inst = InstrumentInfo(imo, "$MOCK_REL/MOCKFT/instrument_info")
		obs = Observation(SpinningScanningStrategy(imo, "$MOCK_REL/scanning_parameters"), inst, dets;
			start_time_s = 100.0, duration_s = 60.0, hwp = IdealHWP(ang_speed_radpsec = inst.hwp_rpm * 2π / 60))
		pt, χ = get_pointings(obs)
		ref = read(j["pointings_ecl"])
		@test size(pt) == size(ref)
		@test maximum(abs.(pt[1, :, :] .- ref[1, :, :])) < 1e-10
		@test maximum(angdiff.(pt[2, :, :], ref[2, :, :])) < 1e-10
		@test maximum(angdiff.(pt[3, :, :], ref[3, :, :])) < 1e-10
		@test maximum(angdiff.(χ, read(j["hwp_angle"]))) < 1e-11
	end
end

@testset "metadata conversion" begin
	det = DetectorInfo(Dict{String, Any}("name" => "x", "pixel" => "12", "squid" => 3.0, "fknee_mhz" => 20,
		"quat" => [0.0, 0.0, 0.0, 2.0], "unknown_key" => 1))
	@test (det.pixel, det.squid, det.fknee_mhz) == (12, 3, 20.0)
	@test det.quat.quats[1] == one(Quaternion{Float64})
	@test_throws ArgumentError DetectorInfo(Dict{String, Any}("pixel" => "A1"))
	@test_throws ArgumentError DetectorInfo(Dict{String, Any}("quat" => Dict("start_time" => "2030-01-01T00:00:00",
		"sampling_rate_hz" => 1.0, "quats" => [[0, 0, 0, 1], [0, 0, 0, 1]])))
	@test_throws ArgumentError DetectorInfo(Dict{String, Any}("quat" => Dict("start_time" => 0.0,
		"sampling_rate_hz" => 1.0, "quaternion_matrix_shape" => [2, 4])))
	ch = FreqChannelInfo(bandcenter_ghz = 140.0)
	@test (ch.channel, ch.number_of_detectors, ch.detector_names) == ("140.0 GHz", 1, ["det0"])
	@test FreqChannelInfo(bandcenter_ghz = 40.0, number_of_detectors = 2).detector_names == ["det0", "det1"]
	@test_throws ArgumentError FreqChannelInfo(bandcenter_ghz = 40.0, number_of_detectors = 2, detector_names = ["a"])
	@test InstrumentInfo(Dict{String, Any}()).name == "mock-instrument"
end

@testset "detector_list_from_parameters" begin
	ch_url = "$MOCK_REL/MOCKFT/M-100/channel_info"
	@test length(detector_list_from_parameters(imo, [Dict("channel_info_obj" => ch_url, "num_of_detectors_from_channel" => 2)])) == 2
	bore = only(detector_list_from_parameters(imo, [Dict("channel_info_obj" => ch_url, "use_only_one_boresight_detector" => true,
		"detector_name" => "bore")]))
	@test bore.name == "bore"
	@test bore.net_ukrts == 50.0
	@test is_constant(bore.quat) && bore.quat.quats[1] == one(Quaternion{Float64})
	over = only(detector_list_from_parameters(imo, [Dict("detector_info_obj" => "$MOCK_REL/MOCKFT/M-100/M100_00_B/detector_info",
		"fwhm_arcmin" => 99.0)]))
	@test (over.name, over.fwhm_arcmin) == ("M100_00_B", 99.0)
	plain = only(detector_list_from_parameters(nothing, [Dict("name" => "p", "sampling_rate_hz" => 5.0)]))
	@test (plain.name, plain.sampling_rate_hz) == ("p", 5.0)
	@test_throws ArgumentError detector_list_from_parameters(nothing, [Dict("channel_info_obj" => ch_url)])
end

# Optional check against a real instrument database: set FALCONS_TEST_IMO to a
# directory containing an (uncompressed) schema.json. The detector quaternions
# are compared with the (θ, φ) listed in pointing_theta_phi_psi_deg.
if haskey(ENV, "FALCONS_TEST_IMO")
	@testset "real IMO ($(ENV["FALCONS_TEST_IMO"]))" begin
		real = Imo(ENV["FALCONS_TEST_IMO"])
		tag = first(sort(collect(keys(real.releases))))
		inst_files = [df for df in values(real.data_files) if real.quantities[df.quantity].name == "instrument_info" && tag in df.release_tags]
		@test !isempty(inst_files)
		inst = InstrumentInfo(real, first(inst_files).uuid)
		ch = FreqChannelInfo(real, first(inst.channel_objs))
		dets = detector_list_from_parameters(real, [Dict("channel_info_obj" => "/data_files/$(first(inst.channel_objs))")])
		@test length(dets) == ch.number_of_detectors
		for det in dets
			(is_constant(det.quat) && det.pointing_theta_phi_psi_deg !== nothing) || continue
			θ, φ, _ = compute_pointing_and_orientation(det.quat.quats[1])
			@test θ ≈ deg2rad(det.pointing_theta_phi_psi_deg[1]) atol = 1e-8
			θ > 1e-6 && @test angdiff(φ, deg2rad(det.pointing_theta_phi_psi_deg[2])) < 1e-8
		end
		scan_files = [df for df in values(real.data_files) if haskey(df.metadata, "spin_sun_angle_deg") && tag in df.release_tags]
		if !isempty(scan_files)
			scan = SpinningScanningStrategy(real, first(scan_files).uuid)
			obs = Observation(scan, inst, dets[1:min(4, end)]; start_time_s = 0.0, duration_s = 10.0)
			pt, _ = get_pointings(obs)
			@test all(isfinite, pt)
		end
	end
end
