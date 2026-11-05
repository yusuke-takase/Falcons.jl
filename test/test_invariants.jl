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
	θ_g, φ_g, _ = Falcons._rotate_coordinates_and_pol_e2g_for_one_sample(0.0, 0.0, 0.0)
	@test rad2deg(φ_g) ≈ 96.3843 atol = 1e-3
	@test rad2deg(π / 2 - θ_g) ≈ 29.8113 atol = 1e-3
	v = rotate_vector_e2g!([0.0, 0.0, 1.0])
	@test v ≈ Falcons.NORTH_POLE_VEC
	@test v ≈ [sin(θ_g) * cos(φ_g), sin(θ_g) * sin(φ_g), cos(θ_g)]
end

@testset "hit maps from angles and pixels" begin
	ss = gen_ScanningStrategy(nside = 8, duration = 300, sampling_rate = 1.0)
	θ, φ, _, _ = get_pointings(ss, 0, 300)
	pix, _, _ = get_pointing_pixels(ss, 0, 300)
	@test hitmap_from_angles(8, θ, φ) == hitmap_from_pixels(8, pix)
	@test sum(hitmap_from_pixels(8, pix)) == 300
end

@testset "known satellites" begin
	@test get_satellite("LiteBIRD").alpha == 45.0
	@test_throws ArgumentError get_satellite("Unknown")
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

@testset "division need not divide the duration" begin
	ss = gen_ScanningStrategy(nside = 8, duration = 601, sampling_rate = 1.0, hwp_rpm = 46.0)
	f1 = get_scanfield(ss; division = 1, spin_n = [0, 2], spin_m = [-4, 4])
	f7 = get_scanfield(ss; division = 7, spin_n = [0, 2], spin_m = [-4, 4])
	@test sum(f7.hitmap) == 601
	@test f7.hitmap == f1.hitmap
	@test isequal(isnan.(f7.h), isnan.(f1.h))
	seen = .!isnan.(f1.h)
	@test maximum(abs.(f7.h[seen] .- f1.h[seen])) < 1e-12
	@test_throws ArgumentError h_nm(f1, 1, 0)
end

@testset "division does not change the samples" begin
	# 0.5 Hz for 10 s: the samples are t = 0, 2, 4, 6, 8 s whatever the chunks.
	ss = gen_ScanningStrategy(nside = 8, duration = 10, sampling_rate = 0.5, hwp_rpm = 46.0)
	f1 = get_scanfield(ss; division = 1, spin_n = [0, 2], spin_m = [-4, 4])
	@test sum(f1.hitmap) == 5
	for division in (2, 3, 5, 7)
		f = get_scanfield(ss; division, spin_n = [0, 2], spin_m = [-4, 4])
		@test f.hitmap == f1.hitmap
		@test isequal(f.h, f1.h)
	end
	f3 = get_scanfield_step(ss; division = 3, spin_n = [0, 2], spin_m = [-4, 4], step = 4)
	@test isequal(f3.h, get_scanfield_step(ss; division = 1, spin_n = [0, 2], spin_m = [-4, 4], step = 4).h)
end
