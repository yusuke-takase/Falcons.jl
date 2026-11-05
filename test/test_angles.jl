@testset "polarization angle convention" begin
	pa = Falcons.polarization_angle
	# At (θ, φ) = (π/2, 0) the local East is +y and the local North is +z.
	@test pa(π / 2, 0.0, [0.0, 1.0, 0.0]) ≈ 0.0 atol = 1e-15
	@test pa(π / 2, 0.0, [0.0, 0.0, 1.0]) ≈ π / 2
	@test pa(π / 2, 0.0, [0.0, 0.0, -1.0]) ≈ -π / 2
	@test abs(pa(π / 2, 0.0, [0.0, -1.0, 0.0])) ≈ π

	# ψ is measured from East (e_φ) towards North (−e_θ), as in litebird_sim.
	θ, φ = 1.1, 2.3
	e_φ = [-sin(φ), cos(φ), 0.0]
	e_θ = [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
	for ψ in range(-3.0, 3.0; length = 13)
		p = cos(ψ) .* e_φ .- sin(ψ) .* e_θ
		@test pa(θ, φ, p) ≈ ψ atol = 1e-12
	end
end

@testset "unit conversions" begin
	@test rpm2angfreq(60.0) ≈ 2π
	@test period2rpm(2.0) == 0.5
	@test period2rpm(30.0; unit = "sec") == 2.0
	@test period2rpm(0.5; unit = "hour") ≈ 1 / 30
	@test rpm2period(0.5) == 2.0
	@test rpm2period(2.0; unit = "sec") == 30.0
	@test arcmin2rad(60.0) ≈ deg2rad(1.0)
end
