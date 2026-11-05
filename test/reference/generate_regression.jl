# Regenerate `regression_master.h5` from the current source tree:
#
#     julia --project=. test/reference/generate_regression.jl
#
# The committed file was produced from master a04df80 (Falcons v1.0.2) before
# any refactoring, so `test_regression.jl` checks that refactors keep the
# numerical output unchanged. Do not regenerate it after behaviour changes.

using Falcons
using HDF5

include(joinpath(@__DIR__, "regression_cases.jl"))

write_compressed(f, name, x) = (f[name, chunk = size(x), deflate = 9] = x)

function main(path = joinpath(@__DIR__, "regression_master.h5"))
	sha = try
		readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
	catch
		"unknown"
	end
	h5open(path, "w") do f
		attributes(f)["generated_from"] = sha
		attributes(f)["julia_version"] = string(VERSION)
		for (name, run) in POINTING_CASES
			θ, φ, ψ, t = run()
			write_compressed(f, "pointings/$name/theta", θ)
			write_compressed(f, "pointings/$name/phi", φ)
			write_compressed(f, "pointings/$name/psi", ψ)
			write_compressed(f, "pointings/$name/time", collect(Float64, t))
		end
		for (name, run) in SCANFIELD_CASES
			field = run()
			write_compressed(f, "scanfield/$name/hitmap", field.hitmap)
			write_compressed(f, "scanfield/$name/h_re", real.(field.h))
			write_compressed(f, "scanfield/$name/h_im", imag.(field.h))
			f["scanfield/$name/spin_n"] = Float64.(field.n)
			f["scanfield/$name/spin_m"] = Float64.(field.m)
		end
	end
	println("wrote ", path)
end

main()
