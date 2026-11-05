@testset "HDF5 scan-field writer" begin
	ss = gen_ScanningStrategy(nside = 8, duration = 120, sampling_rate = 1.0, hwp_rpm = 10.0)
	field = get_scanfield(ss; division = 1, spin_n = [0, 2], spin_m = [-4, 4])
	mktempdir() do dir
		create_h5_file(dir, "/boresight.h5", field)
		h5open(joinpath(dir, "boresight.h5")) do f
			# Keys read by SBM's `ScanFields.load_hdf5`.
			for key in ("hitmap", "h", "toml/spin_n", "toml/spin_m", "ss/nside", "ss/duration", "ss/sampling_rate")
				@test haskey(f, key)
			end
			@test read(f["hitmap"]) == field.hitmap
			@test isequal(read(f["h"]), field.h)
			@test read(f["toml/spin_n"]) == [0.0, 2.0]
			@test read(f["toml/spin_m"]) == [-4.0, 4.0]
			@test read(f["ss/nside"]) == 8
			@test read(f["ss/sampling_rate"]) == 1.0
		end
	end
end
