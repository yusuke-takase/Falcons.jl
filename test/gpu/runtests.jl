# GPU tests, run separately from the main suite because they need CUDA. From
# the repository root (Julia 1.10 or later):
#
#     julia --project=test/gpu -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'   # once
#     julia -t 4 --project=test/gpu test/gpu/runtests.jl
#
# `Pkg.develop` makes the environment use this checkout of Falcons (recorded in
# test/gpu/Manifest.toml, which is not tracked). The GPU results are compared
# with the CPU backend (which is itself checked against litebird_sim in the
# main suite).

import CUDA  # `import`: HDF5 and CUDA both export `attributes`
using Falcons
using HDF5
using Healpix
using Test

realpath(pkgdir(Falcons)) == realpath(joinpath(@__DIR__, "..", "..")) ||
	error("the GPU tests load Falcons from $(pkgdir(Falcons)); run `Pkg.develop(path = \".\")` in the test/gpu environment")

const REFDIR = joinpath(@__DIR__, "..", "reference")
angdiff(a, b) = abs(rem2pi(a - b, RoundNearest))
include(joinpath(@__DIR__, "..", "lbs_reference.jl"))
include(joinpath(@__DIR__, "..", "nonfinite_cases.jl"))

const EXT = Base.get_extension(Falcons, :FalconsCUDAExt)

"Sample times and HWP angles evaluated in a kernel, with the device copy of a plan."
function _k_times!(t, χ, p, n)
	i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
	if i <= n
		@inbounds t[i] = Falcons._plan_time(p, i)
		@inbounds χ[i] = Falcons._hwp_at(p, i)
	end
	return nothing
end

function device_times(obs)
	plan, _ = Falcons._plan_data(obs, Ecliptic, :all)
	dp, _keep = EXT._upload(plan)
	n = obs.n_samples
	t, χ = CUDA.zeros(Float64, n), CUDA.zeros(Float64, n)
	CUDA.@cuda threads = 256 blocks = cld(n, 256) _k_times!(t, χ, dp, n)
	host_t = [Falcons._plan_time(plan, k) for k in 1:n]
	host_χ = [Falcons._hwp_at(plan, k) for k in 1:n]
	return Array(t), Array(χ), host_t, host_χ
end

bits(x) = reinterpret(UInt64, x)

if !CUDA.functional()
	@warn "CUDA is not functional: skipping the GPU tests"
	exit(0)
end

const NGPU = length(CUDA.devices())
const BACKENDS = NGPU > 1 ? [GPUBackend(), GPUBackend(num_gpus = min(NGPU, 4))] : [GPUBackend()]
const ATOL = 1e-12

function compare_pointings(obs)
	for coord in (Ecliptic, Galactic)
		cpu, _ = get_pointings(obs; output_coordinate_system = coord)
		for backend in BACKENDS
			gpu, _ = get_pointings(obs; output_coordinate_system = coord, backend)
			@test size(gpu) == size(cpu)
			@test maximum(abs.(gpu[1, :, :] .- cpu[1, :, :])) < ATOL
			@test maximum(angdiff.(gpu[2, :, :], cpu[2, :, :])) < ATOL
			@test maximum(angdiff.(gpu[3, :, :], cpu[3, :, :])) < ATOL
		end
	end
end

function compare_maps(obs; nside = 32)
	for nest in (false, true)
		cpu = get_pointing_pixels(obs; nside, nest)
		for backend in BACKENDS
			gpu = get_pointing_pixels(obs; nside, nest, backend)
			@test gpu.pixels == cpu.pixels
			@test maximum(angdiff.(gpu.psi, cpu.psi)) < ATOL
			@test get_hitmap(obs; nside, nest, backend) == get_hitmap(obs; nside, nest)
		end
	end
	couples = nm_grid(0:4, (-4, 0, 4))
	cpu = make_h_maps(obs; nside, n_m_couples = couples)
	for backend in BACKENDS
		gpu = make_h_maps(obs; nside, n_m_couples = couples, backend)
		@test gpu.hits == cpu.hits
		@test maximum(abs.(gpu.h .- cpu.h)) < 1e-10
		@test gpu.detectors == cpu.detectors
	end
end

@testset "GPU backend ($(NGPU) devices)" begin
	h5open(LBS_REF) do f
		fs = lbs_cfg(f, "sampling_rate_hz")
		T = lbs_cfg(f, "duration_s")
		t0 = lbs_cfg(f, "start_time_b")

		@testset "analytic attitude, exact and interpolated" begin
			compare_pointings(lbs_observation_b(f))
			compare_pointings(lbs_observation_b(f; delta_time_s = 2.0))
			compare_maps(lbs_observation_b(f))
			compare_maps(lbs_observation_b(f; delta_time_s = 2.0))
			# Also directly against litebird_sim.
			pt, _ = get_pointings(lbs_observation_b(f); backend = GPUBackend())
			dθ, dφ, dψ = lbs_pointing_errors(pt, read(f["B/pointings_ecl"]), read(f["B/sample_index"]))
			@test max(dθ, dφ, dψ) < 1e-10
			@test get_hitmap(lbs_observation_b(f); nside = 32, backend = GPUBackend()) == read(f["F/hitmap"])
		end

		@testset "time-dependent detector quaternion" begin
			g = f["G_exact"]
			det2 = RotQuaternion(read(g["det2_quats"]); start_time_s = read(attributes(g)["det2_start_time_s"]),
				sampling_rate_hz = read(attributes(g)["det2_sampling_rate_hz"]))
			dets = lbs_detectors_b(f)
			dets[3] = DetectorInfo(name = "det2", sampling_rate_hz = fs, quat = det2)
			for delta in (nothing, 2.0)
				obs = Observation(lbs_scan(f), lbs_instrument_b(f), dets; start_time_s = t0, duration_s = T,
					hwp = lbs_hwp_b(f), delta_time_s = delta)
				compare_pointings(obs)
				compare_maps(obs)
			end
		end

		@testset "attitude time series and spin perturbation" begin
			h = f["H"]
			attitude = RotQuaternion(read(h["attitude_quats"]); start_time_s = t0, sampling_rate_hz = 1.0)
			obs = Observation(attitude, lbs_instrument_b(f), lbs_detectors_b(f); start_time_s = t0, duration_s = T,
				hwp = lbs_hwp_b(f))
			compare_pointings(obs)
			compare_maps(obs)
			i = f["I"]
			perturbation = RotQuaternion(read(i["perturbation_quats"]); start_time_s = t0, sampling_rate_hz = fs)
			compare_pointings(lbs_observation_b(f; spin_perturbation = perturbation))
			compare_maps(lbs_observation_b(f; spin_perturbation = perturbation))
		end
	end

	@testset "one day at 19 Hz, split observations" begin
		scan, inst = satellite_preset(:LiteBIRD)
		dets = [DetectorInfo(name = "a", sampling_rate_hz = 19.0), DetectorInfo(name = "b", sampling_rate_hz = 19.0,
			quat = detector_quat(0.05, 1.0; ψ = 0.3))]
		obs = Observation(scan, inst, dets; start_time_s = 0.0, duration_s = 86400.0,
			hwp = IdealHWP(ang_speed_radpsec = 46 * 2π / 60))
		couples = nm_grid(0:4, (-4, 0, 4))
		parts = split_observation(obs, 3)
		cpu = make_h_maps(parts; nside = 128, n_m_couples = couples)
		for backend in BACKENDS
			gpu = make_h_maps(parts; nside = 128, n_m_couples = couples, backend)
			@test gpu.hits == cpu.hits
			@test maximum(abs.(gpu.h .- cpu.h)) < 1e-10
		end
		# Splitting moves the sample times by rounding (~1e-11 s at t ~ 1e4 s),
		# which the m = ±4 HWP terms amplify: a few 1e-10 in h.
		whole = make_h_maps(obs; nside = 128, n_m_couples = couples, backend = GPUBackend())
		@test whole.hits == cpu.hits
		@test maximum(abs.(whole.h .- cpu.h)) < 1e-9
	end

	@testset "times and HWP angles on the device" begin
		# The h-map kernel evaluates the sample times and HWP angles on the
		# device: they must be bit-identical to the host values (explicit fma).
		scan, inst = satellite_preset(:LiteBIRD)
		h5open(LBS_REF) do f
			obs_b = lbs_observation_b(f)
			for obs in (obs_b, split_observation(obs_b, 3)[3], lbs_observation_b(f; delta_time_s = 2.0),
				Observation(scan, inst, DetectorInfo(name = "d", sampling_rate_hz = 19.0); start_time_s = 7.3e7,
					duration_s = 3600.0, start_time_simulation_s = 1.0, hwp = IdealHWP(ang_speed_radpsec = 4.8,
						start_angle_rad = 0.3)))
				t, χ, host_t, host_χ = device_times(obs)
				@test bits.(t) == bits.(host_t)
				@test bits.(χ) == bits.(host_χ)
				@test host_χ == get_pointings(obs)[2]
			end
		end
	end

	@testset "non-finite inputs stop the GPU backend as the CPU one" begin
		for (name, obs) in nonfinite_observations(), backend in BACKENDS
			@test_throws ArgumentError get_pointings(obs; backend)
			@test_throws ArgumentError get_pointing_pixels(obs; nside = 8, backend)
			@test_throws ArgumentError get_hitmap(obs; nside = 8, backend)
			@test_throws ArgumentError make_h_maps(obs; nside = 8, backend)
		end
		# Backstop in the kernels: a plan corrupted after the host checks.
		scan, inst = satellite_preset(:LiteBIRD)
		for delta in (nothing, 2.0)
			obs = Observation(scan, inst, DetectorInfo(name = "a", sampling_rate_hz = 1.0); start_time_s = 0.0,
				duration_s = 10.0, delta_time_s = delta)
			plan, _ = Falcons._plan_data(obs, Galactic, :all)
			plan.dets.constant[1] = Quaternion(NaN, 0.0, 0.0, 1.0)
			resol = Resolution(8)
			npix = resol.numOfPixels
			@test_throws ArgumentError Falcons._cpu_pointing_pixels(plan, obs.n_samples, 1, resol, false)
			for backend in BACKENDS
				@test_throws ArgumentError Falcons._gpu_pointing_pixels(plan, obs.n_samples, 1, resol, false, backend)
				@test_throws ArgumentError Falcons._gpu_hits!(zeros(Int, npix), plan, obs.n_samples, 1, resol, false, backend)
				@test_throws ArgumentError Falcons._gpu_h_maps!(zeros(ComplexF64, 1, npix, 1), zeros(Int, npix, 1), plan,
					obs.n_samples, [1], [2.0], [0.0], resol, backend)
			end
		end
	end

	@testset "the caller's current device is kept" begin
		scan, inst = satellite_preset(:LiteBIRD)
		obs = Observation(scan, inst, DetectorInfo(name = "d", sampling_rate_hz = 1.0); start_time_s = 0.0, duration_s = 100.0)
		for id in (NGPU > 1 ? (NGPU - 1, 0) : (0,))
			CUDA.device!(id)
			for backend in BACKENDS
				make_h_maps(obs; nside = 8, backend)
				get_hitmap(obs; nside = 8, backend)
				get_pointings(obs; backend)
				@test CUDA.deviceid(CUDA.device()) == id
			end
		end
		CUDA.device!(0)
	end

	@testset "errors" begin
		scan, inst = satellite_preset(:LiteBIRD)
		obs = Observation(scan, inst, DetectorInfo(name = "d", sampling_rate_hz = 1.0); start_time_s = 0.0, duration_s = 10.0)
		@test_throws ArgumentError get_pointings(obs; backend = GPUBackend(gpu_ids = [NGPU]))
		short = RotQuaternion([one(Quaternion{Float64}), one(Quaternion{Float64})]; start_time_s = 0.0, sampling_rate_hz = 1.0)
		bad = Observation(short, inst, DetectorInfo(name = "d", sampling_rate_hz = 1.0); start_time_s = 0.0, duration_s = 10.0)
		@test_throws ArgumentError get_pointings(bad; backend = GPUBackend())
	end
end
