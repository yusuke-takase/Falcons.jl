# GPU tests, run separately from the main suite because they need CUDA:
#
#     julia --project=test/gpu -e 'using Pkg; Pkg.instantiate()'   # once
#     julia -t 4 --project=test/gpu test/gpu/runtests.jl
#
# The GPU results are compared with the CPU backend (which is itself checked
# against litebird_sim in the main suite).

import CUDA  # `import`: HDF5 and CUDA both export `attributes`
using Falcons
using HDF5
using Healpix
using Test

const REFDIR = joinpath(@__DIR__, "..", "reference")
angdiff(a, b) = abs(rem2pi(a - b, RoundNearest))
include(joinpath(@__DIR__, "..", "lbs_reference.jl"))

if !CUDA.functional()
	@warn "CUDA is not functional: skipping the GPU tests"
	exit(0)
end

const NGPU = length(CUDA.devices())
const BACKENDS = NGPU > 1 ? [GPUBackend(), GPUBackend(num_gpus = min(NGPU, 4))] : [GPUBackend()]
const ATOL = 1e-12

function compare_pointings(obs)
	for coord in (Ecliptic, Galactic)
		cpu, χ_cpu = get_pointings(obs; output_coordinate_system = coord)
		for backend in BACKENDS
			gpu, χ_gpu = get_pointings(obs; output_coordinate_system = coord, backend)
			@test size(gpu) == size(cpu)
			@test maximum(abs.(gpu[1, :, :] .- cpu[1, :, :])) < ATOL
			@test maximum(angdiff.(gpu[2, :, :], cpu[2, :, :])) < ATOL
			@test maximum(angdiff.(gpu[3, :, :], cpu[3, :, :])) < ATOL
			@test isequal(χ_gpu, χ_cpu)
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

	@testset "errors" begin
		scan, inst = satellite_preset(:LiteBIRD)
		obs = Observation(scan, inst, DetectorInfo(name = "d", sampling_rate_hz = 1.0); start_time_s = 0.0, duration_s = 10.0)
		@test_throws ArgumentError get_pointings(obs; backend = GPUBackend(gpu_ids = [NGPU]))
		short = RotQuaternion([one(Quaternion{Float64}), one(Quaternion{Float64})]; start_time_s = 0.0, sampling_rate_hz = 1.0)
		bad = Observation(short, inst, DetectorInfo(name = "d", sampling_rate_hz = 1.0); start_time_s = 0.0, duration_s = 10.0)
		@test_throws ArgumentError get_pointings(bad; backend = GPUBackend())
	end
end
