@testset "rotation rates and periods" begin
	# 0.05 rpm: one revolution every 20 min = 1200 s = 1/3 h.
	@test rpm2period(0.05) ≈ 20.0
	@test rpm2period(0.05; unit = "sec") ≈ 1200.0
	@test rpm2period(0.05; unit = "hour") ≈ 1 / 3
	for unit in ("min", "sec", "hour")
		@test period2rpm(rpm2period(0.05; unit); unit) ≈ 0.05
	end
	@test period2rpm(192.348) ≈ 1 / 192.348
	@test rpm2angfreq(60.0) ≈ 2π
	@test arcmin2rad(60.0) ≈ deg2rad(1.0)
end
