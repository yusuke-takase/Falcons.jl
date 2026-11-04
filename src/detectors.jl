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

Frequency-channel parameters, with the same fields, defaults and
post-initialisation as `FreqChannelInfo` in litebird_sim (the bandpass object
is not included):

- `channel` defaults to `"<bandcenter_ghz> GHz"` with one decimal;
- `number_of_detectors` is inferred from `detector_names` or `detector_objs`
  when it is 0, and missing names become `det0`, `det1`, …;
- a zero `net_channel_ukrts` or `net_detector_ukrts` is derived from the other
  one with the √N scaling;
- `detector_objs` are stored as UUID strings (a URL keeps its last component).
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

	function FreqChannelInfo(bandcenter_ghz, channel, bandwidth_ghz, band_freqs_ghz, band_weights, net_detector_ukrts,
		net_channel_ukrts, pol_sensitivity_channel_ukarcmin, sampling_rate_hz, fwhm_arcmin, ellipticity, psi_rad,
		fknee_mhz, fmin_hz, alpha, number_of_detectors, detector_names, detector_objs)
		channel = channel === nothing ? @sprintf("%.1f GHz", bandcenter_ghz) : String(channel)
		names = String.(detector_names)
		objs = url_to_uuid.(String.(detector_objs))
		if number_of_detectors > 0
			if isempty(names)
				names = ["det$(i - 1)" for i in 1:number_of_detectors]
			elseif length(names) != number_of_detectors
				throw(ArgumentError("$(length(names)) detector names for $number_of_detectors detectors"))
			end
			isempty(objs) || length(objs) == number_of_detectors ||
				throw(ArgumentError("$(length(objs)) detector objects for $number_of_detectors detectors"))
		elseif !isempty(names)
			number_of_detectors = length(names)
		elseif !isempty(objs)
			number_of_detectors = length(objs)
		else
			number_of_detectors = 1
			names = ["det0"]
		end
		if net_channel_ukrts == 0 && net_detector_ukrts > 0
			net_channel_ukrts = net_detector_ukrts / sqrt(number_of_detectors)
		elseif net_channel_ukrts > 0 && net_detector_ukrts == 0
			net_detector_ukrts = net_channel_ukrts * sqrt(number_of_detectors)
		end
		return new(bandcenter_ghz, channel, bandwidth_ghz, band_freqs_ghz, band_weights, net_detector_ukrts,
			net_channel_ukrts, pol_sensitivity_channel_ukarcmin, sampling_rate_hz, fwhm_arcmin, ellipticity, psi_rad,
			fknee_mhz, fmin_hz, alpha, number_of_detectors, names, objs)
	end
end

"""
    get_boresight_detector(ch::FreqChannelInfo; name = "mock") -> DetectorInfo

A detector looking along the boresight with the averaged parameters of the
channel, as `FreqChannelInfo.get_boresight_detector` in litebird_sim.
"""
function get_boresight_detector(ch::FreqChannelInfo; name::AbstractString = "mock")
	return DetectorInfo(name = name, channel = ch.channel, fwhm_arcmin = ch.fwhm_arcmin, ellipticity = ch.ellipticity,
		psi_rad = ch.psi_rad, net_ukrts = ch.net_detector_ukrts, fknee_mhz = ch.fknee_mhz, fmin_hz = ch.fmin_hz,
		alpha = ch.alpha, sampling_rate_hz = ch.sampling_rate_hz, bandwidth_ghz = ch.bandwidth_ghz,
		bandcenter_ghz = ch.bandcenter_ghz)
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


# ----------------------------------------------------------------------------
# Construction from IMO metadata (the `from_dict` / `from_imo` methods of
# litebird_sim)

_opt_string(x) = x === nothing ? nothing : string(x)
_opt_int(x::Nothing) = nothing
_opt_int(x::Integer) = Int(x)
_opt_int(x::Real) = isinteger(x) ? Int(x) : throw(ArgumentError("expected an integer, got $x"))
function _opt_int(x::AbstractString)
	v = tryparse(Int, strip(x))
	v === nothing && throw(ArgumentError("expected an integer, got \"$x\""))
	return v
end
_opt_floats(x) = x === nothing ? nothing : Float64.(x)

function _quat_from_metadata(q)
	q isa AbstractDict || return RotQuaternion(Float64.(q))
	start_time = get(q, "start_time", nothing)
	start_time isa AbstractString &&
		throw(ArgumentError("time-dependent quaternions with a calendar start time (\"$start_time\") are not supported; use seconds"))
	haskey(q, "quats") ||
		throw(ArgumentError("time-dependent quaternions stored outside the metadata (\"quaternion_matrix_shape\") are not supported"))
	quats = [Quaternion(Float64.(row)) for row in q["quats"]]
	start_time === nothing && return RotQuaternion(only(quats))
	return RotQuaternion(quats; start_time_s = start_time, sampling_rate_hz = q["sampling_rate_hz"])
end

function _detector_field(name::Symbol, value)
	name === :quat && return _quat_from_metadata(value)
	name in (:pixel, :squid) && return _opt_int(value)
	name in (:wafer, :pixtype, :channel, :pol, :orient) && return _opt_string(value)
	name === :name && return string(value)
	name in (:band_freqs_ghz, :band_weights, :pointing_theta_phi_psi_deg, :pointing_u_v) && return _opt_floats(value)
	return Float64(value)
end

"""
    DetectorInfo(metadata::AbstractDict)
    DetectorInfo(imo::Imo, url)

Detector parameters from a dictionary (e.g. the metadata of an IMO
`detector_info` data file), as `DetectorInfo.from_dict` / `from_imo` in
litebird_sim. Keys that are not fields are ignored. `quat` may be a list of 4
numbers or a dictionary with `start_time` (s), `sampling_rate_hz` and `quats`
for a time-dependent rotation. Unlike litebird_sim, `pixel` and `squid` are
converted to integers and no default bandpass is synthesised.
"""
function DetectorInfo(metadata::AbstractDict)
	kw = Dict{Symbol, Any}()
	for f in fieldnames(DetectorInfo)
		key = string(f)
		haskey(metadata, key) && (kw[f] = _detector_field(f, metadata[key]))
	end
	return DetectorInfo(; kw...)
end

DetectorInfo(imo::Imo, url::AbstractString) = DetectorInfo(_metadata(imo, url))

function _with_fields(det::DetectorInfo, overrides::AbstractDict)
	kw = Dict{Symbol, Any}(f => getfield(det, f) for f in fieldnames(DetectorInfo))
	for f in fieldnames(DetectorInfo)
		key = string(f)
		haskey(overrides, key) && (kw[f] = _detector_field(f, overrides[key]))
	end
	return DetectorInfo(; kw...)
end

"""
    FreqChannelInfo(metadata::AbstractDict)
    FreqChannelInfo(imo::Imo, url)

Channel parameters from a dictionary or an IMO `channel_info` data file, as
`FreqChannelInfo.from_dict` / `from_imo` in litebird_sim.
"""
function FreqChannelInfo(metadata::AbstractDict)
	kw = Dict{Symbol, Any}(:bandcenter_ghz => 0.0)
	for f in fieldnames(FreqChannelInfo)
		key = string(f)
		haskey(metadata, key) || continue
		v = metadata[key]
		kw[f] = f === :channel ? _opt_string(v) :
				f === :number_of_detectors ? Int(v) :
				f in (:detector_names, :detector_objs) ? String.(v) :
				f in (:band_freqs_ghz, :band_weights) ? _opt_floats(v) : Float64(v)
	end
	return FreqChannelInfo(; kw...)
end

FreqChannelInfo(imo::Imo, url::AbstractString) = FreqChannelInfo(_metadata(imo, url))

"""
    InstrumentInfo(metadata::AbstractDict)
    InstrumentInfo(imo::Imo, url)

Instrument geometry from a dictionary or an IMO `instrument_info` data file,
as `InstrumentInfo.from_dict` / `from_imo` in litebird_sim: the angles are read
from `boresight_rotangle_deg`, `spin_boresight_angle_deg` and
`spin_rotangle_deg`, the wafers from `wafers` and `waferspace`.
"""
function InstrumentInfo(d::AbstractDict)
	return InstrumentInfo(
		name = string(get(d, "name", "mock-instrument")),
		boresight_rotangle_rad = deg2rad(Float64(get(d, "boresight_rotangle_deg", 0.0))),
		spin_boresight_angle_rad = deg2rad(Float64(get(d, "spin_boresight_angle_deg", 0.0))),
		spin_rotangle_rad = deg2rad(Float64(get(d, "spin_rotangle_deg", 0.0))),
		hwp_rpm = Float64(get(d, "hwp_rpm", 0.0)),
		number_of_channels = Int(get(d, "number_of_channels", 0)),
		channel_names = String.(get(d, "channel_names", String[])),
		channel_objs = url_to_uuid.(String.(get(d, "channel_objs", String[]))),
		wafer_names = String.(get(d, "wafers", String[])),
		wafer_space_cm = Float64(get(d, "waferspace", 0.0)),
	)
end

InstrumentInfo(imo::Imo, url::AbstractString) = InstrumentInfo(_metadata(imo, url))

"""
    detector_list_from_parameters(imo, definitions) -> Vector{DetectorInfo}

Build detectors from a list of dictionaries, e.g. the `[[detectors]]` tables of
a parameter file, as `detector_list_from_parameters` in litebird_sim. Each
entry is one of

- `channel_info_obj = url`: the detectors of a channel, optionally only the
  first `num_of_detectors_from_channel`, or a single boresight detector named
  `detector_name` when `use_only_one_boresight_detector = true`;
- `detector_info_obj = url`: one detector from the IMO;
- neither: a default `DetectorInfo`.

Any `DetectorInfo` field given in an entry (other than with
`channel_info_obj`) overrides the value read from the IMO.
"""
function detector_list_from_parameters(imo::Union{Nothing, Imo}, definitions::AbstractVector)
	result = DetectorInfo[]
	for def in definitions
		def isa AbstractDict || throw(ArgumentError("each detector definition must be a dictionary, got $(typeof(def))"))
		if haskey(def, "channel_info_obj")
			imo === nothing && throw(ArgumentError("channel_info_obj needs an IMO"))
			ch = FreqChannelInfo(imo, def["channel_info_obj"])
			if get(def, "use_only_one_boresight_detector", false)
				push!(result, get_boresight_detector(ch; name = get(def, "detector_name", "mock_boresight")))
				continue
			end
			n = get(def, "num_of_detectors_from_channel", ch.number_of_detectors)
			append!(result, [DetectorInfo(imo, uuid) for uuid in ch.detector_objs[1:n]])
			continue
		end
		det = if haskey(def, "detector_info_obj")
			imo === nothing && throw(ArgumentError("detector_info_obj needs an IMO"))
			DetectorInfo(imo, def["detector_info_obj"])
		else
			DetectorInfo()
		end
		push!(result, _with_fields(det, def))
	end
	return result
end
