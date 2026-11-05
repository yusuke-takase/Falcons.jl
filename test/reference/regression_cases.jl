# Scan configurations whose outputs are pinned by `regression_master.h5`.
#
# The file is shared by `generate_regression.jl` (which writes the reference
# values) and `test_regression.jl` (which compares against them), so both always
# run exactly the same configurations.

using DataFrames

const REG_NSIDE = 16
const REG_DURATION = 600 # s, sampled at 1 Hz

reg_ss(; kw...) = gen_ScanningStrategy(; nside = REG_NSIDE, duration = REG_DURATION, sampling_rate = 1.0, kw...)

_unit(q) = q ./ sqrt(sum(abs2, q))

# Detector quaternions in (x, y, z, w) order, deliberately non-trivial.
const REG_QUATS = [
	_unit([0.10, -0.20, 0.30, 0.90]),
	_unit([-0.05, 0.15, 0.70, 0.60]),
	_unit([0.02, 0.01, -0.40, 0.91]),
]

# Old-style IMO detectors (with a `pixtype` column): one per telescope prefix.
function reg_imo_pixtype_ss(; kw...)
	names = ["000_000_000_QA_040_T", "001_001_001_UB_100_B", "002_002_002_QA_280_T"]
	info = DataFrame(
		name = names,
		pixtype = ["LP1", "MP1", "HP1"],
		orient = ["Q", "U", "45"],
		pol = ["T", "B", "T"],
	)
	return reg_ss(; gamma = 30.0, quat = REG_QUATS, name = names, info = info, kw...)
end

# "Reformation" IMO detectors (name prefix "999", no `pixtype` column).
function reg_imo_999_ss(; kw...)
	names = ["999_000_000_QA_040_T", "999_001_001_UB_100_B"]
	info = DataFrame(name = names, pol_angle_rad = [0.3, 1.2])
	return reg_ss(; gamma = 30.0, quat = REG_QUATS[1:2], name = names, info = info, kw...)
end

function reg_arbitrary_ss(; kw...)
	quats = [ang_to_quat(deg2rad(5.0), deg2rad(30.0)), ang_to_quat(deg2rad(10.0), deg2rad(200.0))]
	return reg_ss(; start_point = "pole", quat = quats, name = ["det_a", "det_b"], kw...)
end

# name => thunk returning (theta, phi, psi, time)
const POINTING_CASES = [
	"boresight_equator" => () -> get_pointings(reg_ss(), 0, REG_DURATION),
	"boresight_pole" => () -> get_pointings(reg_ss(start_point = "pole", start_angle = 0.3), 0, REG_DURATION),
	"boresight_galactic" => () -> get_pointings(reg_ss(coord = "G"), 0, REG_DURATION),
	"imo_pixtype" => () -> get_pointings(reg_imo_pixtype_ss(), 0, REG_DURATION),
	"imo_999" => () -> get_pointings(reg_imo_999_ss(), 0, REG_DURATION),
	"arbitrary_dets" => () -> get_pointings_arbitrary_dets(reg_arbitrary_ss(), 0, REG_DURATION),
	"boresight_step" => () -> get_pointings(reg_ss(start_angle = 0.1), 0, REG_DURATION, 100),
]

# name => thunk returning a scan field
const SCANFIELD_CASES = [
	"boresight_hwp" => () -> get_scanfield(reg_ss(hwp_rpm = 46.0);
		division = 2, spin_n = [0, 1, 2, 4], spin_m = [-4, 0, 4]),
	"imo_pixtype_hwp" => () -> get_scanfield(reg_imo_pixtype_ss(hwp_rpm = 46.0);
		division = 3, spin_n = [0, 2, 4], spin_m = [-4, 4]),
	"imo_999_hwp" => () -> get_scanfield(reg_imo_999_ss(hwp_rpm = 46.0);
		division = 1, spin_n = [0, 2], spin_m = [-4, 4]),
	"galactic_half_spin" => () -> get_scanfield(reg_ss(coord = "G");
		division = 1, spin_n = [-1, 0, 1 / 2, 2], spin_m = [-1 / 2, 1 / 2]),
	"boresight_step_hwp" => () -> get_scanfield_step(reg_ss(hwp_rpm = 46.0);
		division = 2, spin_n = [0, 2], spin_m = [-4, 4], step = 100),
]
