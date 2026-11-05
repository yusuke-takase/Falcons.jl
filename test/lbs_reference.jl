# Builds the Falcons counterparts of the litebird_sim configurations stored in
# reference/lbs_reference.h5 (see reference/generate_lbs_reference.py).
# HDF5 arrays written by h5py appear here with their axes reversed, e.g. the
# (ndet, nsamp, 3) pointings are read as (3, nsamp, ndet).

const LBS_REF = joinpath(REFDIR, "lbs_reference.h5")

lbs_cfg(f, key) = read(attributes(f["config"])[key])

function lbs_scan(f)
	return SpinningScanningStrategy(
		spin_sun_angle_rad = lbs_cfg(f, "spin_sun_angle_rad"),
		precession_rate_hz = lbs_cfg(f, "precession_rate_hz"),
		spin_rate_hz = lbs_cfg(f, "spin_rate_hz"),
	)
end

lbs_instrument_a(f) = InstrumentInfo(name = "instr_a", spin_boresight_angle_rad = lbs_cfg(f, "spin_boresight_angle_rad"))

function lbs_instrument_b(f)
	return InstrumentInfo(
		name = "instr_b",
		spin_boresight_angle_rad = lbs_cfg(f, "spin_boresight_angle_rad"),
		spin_rotangle_rad = lbs_cfg(f, "spin_rotangle_rad_b"),
		boresight_rotangle_rad = lbs_cfg(f, "boresight_rotangle_rad_b"),
	)
end

function lbs_detectors_b(f)
	angles = read(f["config/detector_angles_deg"])
	polang = read(f["config/pol_angles_rad"])
	fs = lbs_cfg(f, "sampling_rate_hz")
	return [DetectorInfo(name = "det$(i - 1)", sampling_rate_hz = fs, pol_angle_rad = polang[i],
		quat = detector_quat(deg2rad(angles[1, i]), deg2rad(angles[2, i]); ψ = deg2rad(angles[3, i]))) for i in axes(angles, 2)]
end

lbs_hwp_b(f) = IdealHWP(ang_speed_radpsec = lbs_cfg(f, "hwp_ang_speed_radpsec"), start_angle_rad = lbs_cfg(f, "hwp_start_angle_rad"))

function lbs_observation_b(f; kw...)
	return Observation(lbs_scan(f), lbs_instrument_b(f), lbs_detectors_b(f);
		start_time_s = lbs_cfg(f, "start_time_b"), duration_s = lbs_cfg(f, "duration_s"), hwp = lbs_hwp_b(f), kw...)
end

"Maximum deviations (θ, φ, ψ) of the stored samples of `pt` from the reference."
function lbs_pointing_errors(pt, ref, sample_index)
	p = pt[:, sample_index .+ 1, :]
	return (maximum(abs.(p[1, :, :] .- ref[1, :, :])), maximum(angdiff.(p[2, :, :], ref[2, :, :])),
		maximum(angdiff.(p[3, :, :], ref[3, :, :])))
end
