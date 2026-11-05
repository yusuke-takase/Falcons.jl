"""
    DetectorInfo(; name = "", sampling_rate_hz = 0.0, quat = identity, pol_angle_rad = 0.0, ...)

Detector parameters, with the same fields and defaults as `DetectorInfo` in
litebird_sim. Only `name`, `sampling_rate_hz`, `quat` and `pol_angle_rad` are
used to compute pointings; the other fields are carried along.

`quat` is the rotation from the detector frame (beam along z, polarisation
reference along x) to the boresight frame. It may be a `Quaternion`, a
4-element vector `(x, y, z, w)` or a [`RotQuaternion`](@ref) (constant or
time-dependent); it is always stored as a `RotQuaternion`. See
[`detector_quat`](@ref) to place a detector on the focal plane.

`pol_angle_rad` is the polarisation angle relative to the detector x axis. It
is not added to the orientation angle ψ returned by [`get_pointings`](@ref) nor
used in [`make_h_maps`](@ref), as in litebird_sim.
"""
Base.@kwdef struct DetectorInfo
	name::String = ""
	wafer::Union{Nothing, String} = nothing
	pixel::Union{Nothing, Int} = nothing
	pixtype::Union{Nothing, String} = nothing
	channel::Union{Nothing, String} = nothing
	squid::Union{Nothing, Int} = nothing
	sampling_rate_hz::Float64 = 0.0
	fwhm_arcmin::Float64 = 0.0
	ellipticity::Float64 = 1.0
	psi_rad::Float64 = 0.0
	bandcenter_ghz::Float64 = 0.0
	bandwidth_ghz::Float64 = 0.0
	band_freqs_ghz::Union{Nothing, Vector{Float64}} = nothing
	band_weights::Union{Nothing, Vector{Float64}} = nothing
	net_ukrts::Float64 = 0.0
	pol_sensitivity_ukarcmin::Float64 = 0.0
	fknee_mhz::Float64 = 0.0
	fmin_hz::Float64 = 0.0
	alpha::Float64 = 0.0
	pol::Union{Nothing, String} = nothing
	orient::Union{Nothing, String} = nothing
	quat::RotQuaternion{Float64} = RotQuaternion(one(Quaternion{Float64}))
	pol_angle_rad::Float64 = 0.0
	pol_efficiency::Float64 = 1.0
	pointing_theta_phi_psi_deg::Union{Nothing, Vector{Float64}} = nothing
	pointing_u_v::Union{Nothing, Vector{Float64}} = nothing
	g_one_over_k::Float64 = 0.0
	amplitude_2f_k::Float64 = 0.0

	function DetectorInfo(name, wafer, pixel, pixtype, channel, squid, sampling_rate_hz, fwhm_arcmin, ellipticity,
		psi_rad, bandcenter_ghz, bandwidth_ghz, band_freqs_ghz, band_weights, net_ukrts, pol_sensitivity_ukarcmin,
		fknee_mhz, fmin_hz, alpha, pol, orient, quat, pol_angle_rad, pol_efficiency, pointing_theta_phi_psi_deg,
		pointing_u_v, g_one_over_k, amplitude_2f_k)
		q = RotQuaternion(quat)
		q = q isa RotQuaternion{Float64} ? q : _as_float64(q)
		return new(name, wafer, pixel, pixtype, channel, squid, sampling_rate_hz, fwhm_arcmin, ellipticity, psi_rad,
			bandcenter_ghz, bandwidth_ghz, band_freqs_ghz, band_weights, net_ukrts, pol_sensitivity_ukarcmin, fknee_mhz,
			fmin_hz, alpha, pol, orient, q, pol_angle_rad, pol_efficiency, pointing_theta_phi_psi_deg, pointing_u_v,
			g_one_over_k, amplitude_2f_k)
	end
end

function _as_float64(r::RotQuaternion)
	quats = [Quaternion(Float64(q.x), Float64(q.y), Float64(q.z), Float64(q.w)) for q in r.quats]
	is_constant(r) && return RotQuaternion(quats[1])
	return RotQuaternion(quats; start_time_s = r.start_time_s, sampling_rate_hz = r.sampling_rate_hz)
end

"""
    detector_quat(θ, φ; ψ = 0.0) -> Quaternion

Detector-to-boresight rotation `Rz(φ) · Ry(θ) · Rz(ψ)` of a detector whose beam
points at the focal-plane position `(θ, φ)` (rad, `θ` measured from the
boresight) with its polarisation reference rotated by `ψ` around the beam.
"""
detector_quat(θ, φ; ψ = 0.0) = quat_rotation_z(φ) * (quat_rotation_y(θ) * quat_rotation_z(ψ))

"""
    FreqChannelInfo(; bandcenter_ghz, channel = nothing, ...)

Frequency-channel parameters, with the same fields and defaults as
`FreqChannelInfo` in litebird_sim (the bandpass object is not included).
"""
Base.@kwdef struct FreqChannelInfo
	bandcenter_ghz::Float64
	channel::Union{Nothing, String} = nothing
	bandwidth_ghz::Float64 = 0.0
	band_freqs_ghz::Union{Nothing, Vector{Float64}} = nothing
	band_weights::Union{Nothing, Vector{Float64}} = nothing
	net_detector_ukrts::Float64 = 0.0
	net_channel_ukrts::Float64 = 0.0
	pol_sensitivity_channel_ukarcmin::Float64 = 0.0
	sampling_rate_hz::Float64 = 0.0
	fwhm_arcmin::Float64 = 0.0
	ellipticity::Float64 = 1.0
	psi_rad::Float64 = 0.0
	fknee_mhz::Float64 = 0.0
	fmin_hz::Float64 = 1e-5
	alpha::Float64 = 1.0
	number_of_detectors::Int = 0
	detector_names::Vector{String} = String[]
	detector_objs::Vector{String} = String[]
end

"""
    InstrumentInfo(; name = "", boresight_rotangle_rad = 0.0, spin_boresight_angle_rad = 0.0,
                   spin_rotangle_rad = 0.0, hwp_rpm = 0.0, ...)

Instrument geometry, with the same fields and defaults as `InstrumentInfo` in
litebird_sim. The boresight-to-spin rotation is

    bore2spin_quat = Rz(spin_rotangle_rad) · Ry(spin_boresight_angle_rad) · Rz(boresight_rotangle_rad)

`hwp_rpm` is informative only: pass an [`IdealHWP`](@ref) to the
[`Observation`](@ref) to model a rotating half-wave plate.
"""
Base.@kwdef struct InstrumentInfo
	name::String = ""
	boresight_rotangle_rad::Float64 = 0.0
	spin_boresight_angle_rad::Float64 = 0.0
	spin_rotangle_rad::Float64 = 0.0
	hwp_rpm::Float64 = 0.0
	number_of_channels::Int = 0
	channel_names::Vector{String} = String[]
	channel_objs::Vector{String} = String[]
	wafer_names::Vector{String} = String[]
	wafer_space_cm::Float64 = 0.0
end

"""
    bore2spin_quat(inst::InstrumentInfo) -> Quaternion

Rotation from the boresight frame to the spin-axis frame.
"""
function bore2spin_quat(inst::InstrumentInfo)
	q = quat_rotation_z(inst.boresight_rotangle_rad)
	q = quat_rotation_y(inst.spin_boresight_angle_rad) * q
	q = quat_rotation_z(inst.spin_rotangle_rad) * q
	return normalize(q)
end
