# Inputs that the pointing functions must reject, shared by the CPU and GPU tests.

"Observations whose inputs contain a non-finite value or a zero quaternion."
function nonfinite_observations()
	scan, inst = satellite_preset(:LiteBIRD)
	good = DetectorInfo(name = "a", sampling_rate_hz = 1.0)
	gap = [quat_rotation_z(0.01 * k) for k in 0:20]
	gap[5] = Quaternion(NaN, NaN, NaN, NaN)  # a missing sample
	series = RotQuaternion(gap; start_time_s = 0.0, sampling_rate_hz = 1.0)
	obs(attitude = scan, instrument = inst, dets = good; kw...) =
		Observation(attitude, instrument, dets; start_time_s = 0.0, duration_s = 10.0, kw...)
	return [
		"zero detector quaternion" => obs(scan, inst, [good, DetectorInfo(name = "z", sampling_rate_hz = 1.0, quat = [0, 0, 0, 0])]),
		"missing attitude sample" => obs(series),
		"missing perturbation sample" => obs(; spin_perturbation = series),
		"NaN scan parameter" => obs(SpinningScanningStrategy(spin_sun_angle_rad = 0.7, precession_rate_hz = NaN, spin_rate_hz = 1e-3)),
		"NaN instrument angle" => obs(scan, InstrumentInfo(spin_boresight_angle_rad = NaN)),
		"NaN start time" => Observation(scan, inst, good; start_time_s = NaN, duration_s = 10.0),
		"NaN HWP speed" => obs(; hwp = IdealHWP(ang_speed_radpsec = NaN)),
		"infinite delta_time_s" => obs(; delta_time_s = Inf),
	]
end
