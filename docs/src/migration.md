# Migrating from Falcons 1.0

This version follows litebird_sim's conventions and names, which changes both the API and some results:

- ψ and the HWP angle follow litebird_sim (see [Conventions](@ref)); the detector polarisation angle is
  no longer added to the HWP angle.
- h-maps use `exp(+i(nψ + mχ))`, the complex conjugate of the former `exp(−i(nψ + mϕ))`, and unobserved
  pixels are 0 instead of `NaN`. The `(0, 0)` couple holds the hit counts.
- The focal-plane rotation is `bore2spin = Rz(spin_rotangle) · Ry(β) · Rz(boresight_rotangle)`; the
  former `gamma`, `start_point` and name-dependent IMO rules are gone.
- The Ecliptic-to-Galactic rotation uses litebird_sim's values (they differ from the former ones by up to
  1e-7).

## Names

| Falcons 1.0 | Now |
|---|---|
| `gen_ScanningStrategy(...)`, `ScanningStrategy` | [`SpinningScanningStrategy`](@ref), [`InstrumentInfo`](@ref), [`DetectorInfo`](@ref), [`IdealHWP`](@ref) and [`Observation`](@ref) |
| `ss.alpha`, `ss.beta` (deg) | `spin_sun_angle_rad`, `InstrumentInfo.spin_boresight_angle_rad` |
| `ss.prec_rpm`, `ss.spin_rpm` | `precession_rate_hz`, `spin_rate_hz` (rpm / 60) |
| `ss.hwp_rpm` | `IdealHWP(ang_speed_radpsec = 2π · rpm / 60)` |
| `ss.start_point = "pole"` | `InstrumentInfo(spin_rotangle_rad = π)` |
| `ss.start_angle` | `start_time_s` (the anti-Sun direction is at `Ω · t`) |
| `ss.coord = "G"` | `output_coordinate_system = Galactic` |
| `ss.duration`, `ss.sampling_rate` | `Observation(...; duration_s)`, `DetectorInfo.sampling_rate_hz` |
| `ss.quat`, `ss.name`, `ss.info` | `DetectorInfo(name, quat, ...)` |
| `get_satellite("LiteBIRD")` | [`satellite_preset`](@ref)`(:LiteBIRD)` |
| `show_ss(ss)` | `display(scan)` |
| `get_pointings(ss, start, stop)` | [`get_pointings`](@ref)`(obs)`; split long scans with [`split_observation`](@ref) |
| `get_pointings(ss, start, stop, step)`, `get_scanfield_step` | removed |
| `get_pointings_arbitrary_dets`, `ang_to_quat` | `DetectorInfo(quat = detector_quat(θ, φ))` |
| `get_pointing_pixels(ss, start, stop)` | [`get_pointing_pixels`](@ref)`(obs; nside)` |
| `angtod2hitmap`, `pixtod2hitmap` | [`get_hitmap`](@ref) |
| `get_scanfield`, `scanfield`, `h_nm`, `get_hnm_quantify` | [`make_h_maps`](@ref), [`HMapsResult`](@ref), [`h_map`](@ref), [`hitmap`](@ref) |
| `create_h5_file` | [`save_h_maps`](@ref) (litebird_sim layout), [`load_h_maps_from_file`](@ref) |
| `sim_det_scanfields(toml)` | [`run_h_maps_pipeline`](@ref) (litebird_sim parameter names) |
| `gen_imo`, `imo_name!`, `imo_channel!`, `imo_telescope!`, `get_channel_info`, `get_instrument_info`, `get_pol_angle` | [`Imo`](@ref), [`query`](@ref), `DetectorInfo(imo, url)`, `FreqChannelInfo(imo, url)`, `InstrumentInfo(imo, url)`, [`detector_list_from_parameters`](@ref) |
| `rotate_quat` | [`quat_rotation`](@ref) |
| `rotate_vec_ecl_to_gal_coordinates!` | `rotate_vector_e2g!` |
| `binned_mapmake`, `sim_pointing_systematics`, `get_psi_database` and the offset types | removed |
