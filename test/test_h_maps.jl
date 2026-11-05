const LBS_ATOL_H = 1e-10

h5open(LBS_REF) do f
	obs = lbs_observation_b(f)
	for (case, coord) in (("E_gal", Galactic), ("E_ecl", Ecliptic))
		@testset "$case: h-maps against litebird_sim" begin
			g = f[case]
			couples = read(g["n_m_couples"])  # (2, ncouples)
			res = make_h_maps(obs; nside = read(attributes(g)["nside"]), n_m_couples = permutedims(couples),
				output_coordinate_system = coord)
			@test res.coordinate_system == coord
			for det in res.detectors, c in axes(couples, 2)
				n, m = couples[1, c], couples[2, c]
				ref = complex.(read(g["$det/$n,$m/Re"]), read(g["$det/$n,$m/Im"]))
				hm = h_map(res, det, n, m)
				if (n, m) == (0, 0)
					@test hm == ref
					@test hitmap(res, det) == round.(Int, real.(ref))
				else
					@test maximum(abs.(hm .- ref)) < LBS_ATOL_H
				end
			end
		end
	end
end

scan_h, inst_h = satellite_preset(:LiteBIRD)
dets_h = [DetectorInfo(name = "a", sampling_rate_hz = 4.0), DetectorInfo(name = "b", sampling_rate_hz = 4.0, quat = detector_quat(0.1, 0.2; ψ = 0.3))]
obs_h = Observation(scan_h, inst_h, dets_h; start_time_s = 0.0, duration_s = 1800.0, hwp = IdealHWP(ang_speed_radpsec = 2.0))
couples_h = nm_grid(0:2, (-4, 0, 4))

@testset "definition" begin
	res = make_h_maps(obs_h; nside = 16, n_m_couples = couples_h)
	px = get_pointing_pixels(obs_h; nside = 16)
	@test size(res.h) == (length(couples_h), nside2npix(16), 2)
	@test res.duration_s == 1800.0
	@test res.sampling_rate_hz == 4.0
	@test (res.detector_split, res.time_split) == ("full", "full")
	for (d, det) in enumerate(res.detectors)
		hits = hitmap(res, det)
		@test sum(hits) == obs_h.n_samples
		p = px.pixels[1, d]
		in_p = findall(==(p), px.pixels[:, d])
		for (n, m) in couples_h
			expected = (n, m) == (0, 0) ? complex(length(in_p)) :
					   sum(cis.(n .* px.psi[in_p, d] .+ m .* px.hwp_angle[in_p])) / length(in_p)
			@test h_map(res, det, n, m)[p] ≈ expected
		end
		@test all(iszero, h_map(res, det, 2, 4)[hits .== 0])
	end
	@test_throws ArgumentError h_map(res, "a", 3, 0)
	@test_throws ArgumentError h_map(res, "zz", 0, 0)
end

@testset "observations are accumulated" begin
	whole = make_h_maps(obs_h; nside = 16, n_m_couples = couples_h)
	parts = make_h_maps(split_observation(obs_h, 4); nside = 16, n_m_couples = couples_h)
	@test parts.hits == whole.hits
	@test maximum(abs.(parts.h .- whole.h)) < 1e-12
	@test parts.duration_s ≈ whole.duration_s
	# Observations with different detectors are matched by name.
	only_b = Observation(scan_h, inst_h, dets_h[2:2]; start_time_s = 1800.0, duration_s = 600.0)
	mixed = make_h_maps([obs_h, only_b]; nside = 16, n_m_couples = [(0, 0)])
	@test sum(hitmap(mixed, "a")) == obs_h.n_samples
	@test sum(hitmap(mixed, "b")) == obs_h.n_samples + only_b.n_samples
end

@testset "couples" begin
	@test nm_grid([0, 2], [-4, 4]) == [(0, -4), (0, 4), (2, -4), (2, 4)]
	@test make_h_maps(obs_h; nside = 8, n_m_couples = [0 0; 2 0]).n_m_couples == [(0, 0), (2, 0)]
	@test_throws ArgumentError make_h_maps(obs_h; nside = 8, n_m_couples = [(0.5, 0)])
	@test_throws ArgumentError make_h_maps(obs_h; nside = 8, n_m_couples = [(2, 0), (2, 0)])
	res = make_h_maps(obs_h; nside = 8, n_m_couples = [(2, 0)])
	@test sum(hitmap(res, "a")) == obs_h.n_samples  # hit counts are kept without (0, 0)
end

@testset "save and load (litebird_sim layout)" begin
	res = make_h_maps(obs_h; nside = 8, n_m_couples = [(0, 0), (2, 0), (2, -4)])
	mktempdir() do dir
		paths = save_h_maps(res, dir)
		@test basename.(paths) == ["h_maps_det_a.h5", "h_maps_det_b.h5"]
		h5open(paths[1]) do f
			attrs = attributes(f)
			@test read(attrs["coordinate_system"]) == "Galactic"
			@test read(attrs["det"]) == "a"
			@test read(attrs["duration_s"]) == 1800.0
			@test read(attrs["sampling_rate_Hz"]) == 4.0
			@test (read(attrs["detector_split"]), read(attrs["time_split"])) == ("full", "full")
			@test sort(keys(f)) == ["0,0", "2,-4", "2,0"]
			@test read(f["2,-4/Re"]) == real.(h_map(res, "a", 2, -4))
			# Scalar variable-length UTF-8 strings, which h5py (and so
			# litebird_sim) reads as `str` rather than `bytes`.
			for key in ("coordinate_system", "det", "detector_split", "time_split")
				attr = attrs[key]
				dtype = HDF5.datatype(attr)
				@test HDF5.API.h5t_get_class(dtype) == HDF5.API.H5T_STRING
				@test HDF5.API.h5t_is_variable_str(dtype)
				@test HDF5.API.h5t_get_cset(dtype) == HDF5.API.H5T_CSET_UTF8
				@test HDF5.API.h5s_get_simple_extent_type(HDF5.dataspace(attr)) == HDF5.API.H5S_SCALAR
				close(dtype)
			end
		end
		back = load_h_maps_from_file(paths[2])
		@test back.detectors == ["b"]
		@test back.nside == 8
		@test back.coordinate_system == Galactic
		for (n, m) in res.n_m_couples
			@test h_map(back, "b", n, m) == h_map(res, "b", n, m)
		end
		@test hitmap(back, "b") == hitmap(res, "b")
	end
end
