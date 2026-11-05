# The device-friendly pixelisation must agree exactly with Healpix.jl.

@testset "ang2pix ports agree with Healpix.jl" begin
	edges_θ = [0.0, 1e-12, 0.005, 0.01, acos(2 / 3), acos(2 / 3) + 1e-15, π / 2, acos(-2 / 3), π - 0.01, π - 1e-12, π]
	edges_φ = [0.0, 1e-15, π / 4, π / 2, π, 3π / 2, prevfloat(2π)]
	rng_θ = acos.(range(-1, 1; length = 997))
	rng_φ = range(0, prevfloat(2π); length = 1009)
	for nside in (1, 2, 4, 64, 1024)
		resol = Resolution(nside)
		for (θs, φs) in ((edges_θ, edges_φ), (rng_θ, rng_φ))
			@test all(Falcons._ang2pix_ring(resol, θ, φ) == ang2pixRing(resol, θ, φ) for θ in θs, φ in φs)
			@test all(Falcons._ang2pix_nest(resol, θ, φ) == ang2pixNest(resol, θ, φ) for θ in θs, φ in φs)
		end
	end
	@test all(Falcons._spread_bits(v) == Healpix.spreadbits(v) for v in (0, 1, 2, 3, 255, 256, 12345, 2^20 - 1, 2^29 + 7))
end

@testset "_rem2pi" begin
	for x in Float64[0.0, 1.0, -3.0, π, -π, 10.0, 1e3, 12345.678, -98765.4321, 1.5e8, 3e9]
		r = Falcons._rem2pi(x)
		@test abs(r) <= π + 1e-15
		@test angdiff(r, rem2pi(x, RoundNearest)) < 1e-15 * max(1.0, abs(x) / 1e6)
	end
end

@testset "backends" begin
	@test GPUBackend().num_gpus == 1
	@test GPUBackend(num_gpus = 3).gpu_ids === nothing
	@test GPUBackend(gpu_ids = [2, 0]).num_gpus == 2
	@test_throws ArgumentError GPUBackend(num_gpus = 0)
	@test_throws ArgumentError GPUBackend(gpu_ids = [1, 1])
	if Base.get_extension(Falcons, :FalconsCUDAExt) === nothing
		scan, inst = satellite_preset(:LiteBIRD)
		obs = Observation(scan, inst, DetectorInfo(name = "d", sampling_rate_hz = 1.0); start_time_s = 0.0, duration_s = 10.0)
		@test_throws ArgumentError get_pointings(obs; backend = GPUBackend())
	end
end
