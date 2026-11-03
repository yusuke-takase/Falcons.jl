# Checks against values that follow analytically from the scan geometry.

@testset "boresight direction at t = 0" begin
	for (α, β) in ((45.0, 50.0), (30.0, 65.0), (65.0, 30.0))
		ss = gen_ScanningStrategy(nside = 16, duration = 10, sampling_rate = 1.0, alpha = α, beta = β)
		θ, φ, ψ, _ = get_pointings(ss, 0, 10)
		@test θ[1, 1] ≈ deg2rad(90 - α + β)
		@test angdiff(φ[1, 1], 0.0) < 1e-12
		@test ψ[1, 1] ≈ -π / 2
	end
end

@testset "angle ranges" begin
	ss = gen_ScanningStrategy(nside = 16, duration = 3600, sampling_rate = 1.0)
	θ, φ, ψ, _ = get_pointings(ss, 0, 3600)
	@test all(0 .<= θ .<= π)
	@test all(0 .<= φ .< 2π)
	@test all(-π .< ψ .<= π)
end

@testset "one spin period with precession off" begin
	# With no precession, Q(t) = Rz(λ0 + Ω t) R_spin(ω t); after one spin period
	# the pointing equals the t = 0 pointing of a scan started at λ0 + Ω T.
	T = 1200 # s, spin period for 0.05 rpm
	Ω = 2π / (60 * 60 * 24 * 365.25)
	ss1 = gen_ScanningStrategy(nside = 16, duration = T + 1, sampling_rate = 1.0, prec_rpm = 0.0)
	ss2 = gen_ScanningStrategy(nside = 16, duration = 1, sampling_rate = 1.0, prec_rpm = 0.0, start_angle = Ω * T)
	θ1, φ1, ψ1, t1 = get_pointings(ss1, 0, T + 1)
	θ2, φ2, ψ2, _ = get_pointings(ss2, 0, 1)
	@test t1[end] == T
	@test θ1[end, 1] ≈ θ2[1, 1] atol = 1e-10
	@test angdiff(φ1[end, 1], φ2[1, 1]) < 1e-10
	@test angdiff(ψ1[end, 1], ψ2[1, 1]) < 1e-10
end

@testset "ecliptic to galactic" begin
	# The north ecliptic pole sits at (l, b) ≈ (96.3843°, 29.8113°).
	θ_g, φ_g, _ = _rotate_coordinates_and_pol_e2g_for_one_sample(0.0, 0.0, 0.0)
	@test rad2deg(φ_g) ≈ 96.3843 atol = 1e-3
	@test rad2deg(π / 2 - θ_g) ≈ 29.8113 atol = 1e-3
end

@testset "scan field invariants" begin
	ss = gen_ScanningStrategy(nside = 8, duration = 600, sampling_rate = 1.0, hwp_rpm = 46.0)
	spin_n, spin_m = [0, 1, 2], [-4, 0, 4]
	field = get_scanfield(ss; division = 2, spin_n = spin_n, spin_m = spin_m)
	@test sum(field.hitmap) == 600
	seen = field.hitmap .> 0
	@test all(h_nm(field, 0, 0)[seen] .≈ 1)
	@test all(isnan, h_nm(field, 0, 0)[.!seen])
	for n in spin_n, m in spin_m
		@test all(abs.(h_nm(field, n, m)[seen]) .<= 1 + 1e-12)
	end
end
