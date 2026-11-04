const _DURATION_UNITS = (
	("years", 365.25 * 86400.0), ("year", 365.25 * 86400.0), ("days", 86400.0), ("day", 86400.0),
	("hours", 3600.0), ("hour", 3600.0), ("minutes", 60.0), ("min", 60.0), ("sec", 1.0), ("s", 1.0),
)

# A duration in seconds, or a string such as "1 year", "2 days", "30 min" (as litebird_sim).
_duration_s(x::Real) = Float64(x)
function _duration_s(x::AbstractString)
	m = match(r"^\s*([0-9.eE+-]+)\s*([a-z]+)\s*$", x)
	if m !== nothing
		for (unit, factor) in _DURATION_UNITS
			m.captures[2] == unit && return parse(Float64, m.captures[1]) * factor
		end
	end
	throw(ArgumentError("cannot parse the duration \"$x\"; use seconds or e.g. \"1 year\", \"2 days\", \"30 min\""))
end

_start_time_s(x::Real) = Float64(x)
_start_time_s(x) = throw(ArgumentError("start_time must be a number of seconds, got $(repr(x))"))

function _coordinate_system(name::AbstractString)
	name == "Galactic" && return Galactic
	name == "Ecliptic" && return Ecliptic
	throw(ArgumentError("unknown coordinate system \"$name\"; use \"Galactic\" or \"Ecliptic\""))
end

function _section(params, name)
	section = get(params, name, Dict{String, Any}())
	section isa AbstractDict || throw(ArgumentError("[$name] must be a table"))
	return section
end

function _needs_imo(imo, section)
	imo === nothing && throw(ArgumentError("[$section] uses imo_url, but no [imo] flatfile_location is given"))
	return imo
end

function _scanning_strategy_from_parameters(params, imo)
	s = _section(params, "scanning_strategy")
	haskey(s, "imo_url") && return SpinningScanningStrategy(_needs_imo(imo, "scanning_strategy"), s["imo_url"])
	return SpinningScanningStrategy(spin_sun_angle_rad = s["spin_sun_angle_rad"], precession_rate_hz = s["precession_rate_hz"],
		spin_rate_hz = s["spin_rate_hz"])
end

function _instrument_from_parameters(params, imo)
	s = _section(params, "instrument")
	haskey(s, "imo_url") && return InstrumentInfo(_needs_imo(imo, "instrument"), s["imo_url"])
	kw = Dict{Symbol, Any}()
	for f in (:name, :boresight_rotangle_rad, :spin_boresight_angle_rad, :spin_rotangle_rad, :hwp_rpm)
		haskey(s, string(f)) && (kw[f] = s[string(f)])
	end
	return InstrumentInfo(; kw...)
end

function _hwp_from_parameters(params)
	haskey(params, "hwp") || return nothing
	s = _section(params, "hwp")
	return IdealHWP(ang_speed_radpsec = s["ang_speed_radpsec"], start_angle_rad = get(s, "start_angle_rad", 0.0))
end

"""
    run_h_maps_pipeline(parameter_file) -> HMapsResult

Compute and save h-maps from a TOML parameter file whose keys follow
litebird_sim's names:

```toml
[simulation]
base_path = "output"           # output root
start_time = 0.0               # s
duration_s = "1 year"          # seconds, or "N years/days/hours/min/s"

[imo]                          # optional
flatfile_location = "/path/to/imo"

[scanning_strategy]
imo_url = "/releases/<tag>/…/scanning_parameters"
# or spin_sun_angle_rad, precession_rate_hz, spin_rate_hz
# delta_time_s = 2.0           # optional: interpolate as litebird_sim

[instrument]
imo_url = "/releases/<tag>/…/instrument_info"
# or boresight_rotangle_rad, spin_boresight_angle_rad, spin_rotangle_rad

[hwp]                          # optional
ang_speed_radpsec = 4.817      # 2π · rpm / 60
start_angle_rad = 0.0

[[detectors]]                  # see detector_list_from_parameters
channel_info_obj = "/releases/<tag>/…/channel_info"

[h_maps]
nside = 128
n_m_couples = [[0, 0], [2, 0], [4, 0]]
output_coordinate_system = "Galactic"
output_directory = "h_maps"    # relative to base_path
```

The maps are written with [`save_h_maps`](@ref), one file per detector.
"""
function run_h_maps_pipeline(parameter_file::AbstractString)
	params = TOML.parsefile(parameter_file)
	sim = _section(params, "simulation")
	haskey(sim, "duration_s") || throw(ArgumentError("[simulation] needs duration_s"))
	imo = haskey(params, "imo") ? Imo(_section(params, "imo")["flatfile_location"]) : nothing

	scan = _scanning_strategy_from_parameters(params, imo)
	inst = _instrument_from_parameters(params, imo)
	dets = detector_list_from_parameters(imo, get(params, "detectors", []))
	isempty(dets) && throw(ArgumentError("no [[detectors]] given"))
	delta_time_s = get(_section(params, "scanning_strategy"), "delta_time_s", nothing)
	obs = Observation(scan, inst, dets; start_time_s = _start_time_s(get(sim, "start_time", 0.0)),
		duration_s = _duration_s(sim["duration_s"]), hwp = _hwp_from_parameters(params), delta_time_s)

	h = _section(params, "h_maps")
	haskey(h, "nside") || throw(ArgumentError("[h_maps] needs nside"))
	couples = [(c[1], c[2]) for c in get(h, "n_m_couples", [[0, 0], [2, 0], [4, 0]])]
	res = make_h_maps(obs; nside = h["nside"], n_m_couples = couples,
		output_coordinate_system = _coordinate_system(get(h, "output_coordinate_system", "Galactic")))
	save_h_maps(res, joinpath(get(sim, "base_path", "."), get(h, "output_directory", "h_maps")))
	return res
end
