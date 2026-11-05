# Regression against the outputs of master a04df80 (see
# `reference/generate_regression.jl`). These tests guard refactors that must
# not change any number; they are replaced by litebird_sim reference values
# once the conventions change.

include(joinpath(REFDIR, "regression_cases.jl"))

const REG_ATOL_ANGLE = 1e-10 # rad
const REG_ATOL_H = 1e-10

@testset "regression vs master a04df80" begin
	h5open(joinpath(REFDIR, "regression_master.h5")) do f
		for (name, run) in POINTING_CASES
			@testset "pointings/$name" begin
				θ, φ, ψ, t = run()
				g = f["pointings/$name"]
				θ_ref, φ_ref, ψ_ref = read(g["theta"]), read(g["phi"]), read(g["psi"])
				@test size(θ) == size(θ_ref)
				@test maximum(abs.(θ .- θ_ref)) < REG_ATOL_ANGLE
				@test maximum(angdiff.(φ, φ_ref)) < REG_ATOL_ANGLE
				@test maximum(angdiff.(ψ, ψ_ref)) < REG_ATOL_ANGLE
				@test collect(Float64, t) ≈ read(g["time"])
			end
		end
		for (name, run) in SCANFIELD_CASES
			@testset "scanfield/$name" begin
				field = run()
				g = f["scanfield/$name"]
				@test Float64.(field.n) == read(g["spin_n"])
				@test Float64.(field.m) == read(g["spin_m"])
				@test field.hitmap == read(g["hitmap"])
				h_ref = complex.(read(g["h_re"]), read(g["h_im"]))
				@test size(field.h) == size(h_ref)
				@test isnan.(field.h) == isnan.(h_ref)
				seen = .!isnan.(h_ref)
				@test maximum(abs.(field.h[seen] .- h_ref[seen])) < REG_ATOL_H
			end
		end
	end
end
