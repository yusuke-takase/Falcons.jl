# Instrument database and parameter files

## Reading an IMO

Falcons reads instrument models (IMO) in the flat-file layout used by litebird_sim (via libinsdb): a
`schema.json` with entities, quantities, data files and releases.

```julia
using Falcons

imo = Imo("/path/to/imo")                               # directory containing schema.json
rel = "/releases/<tag>/satellite"

scan = SpinningScanningStrategy(imo, "$rel/scanning_parameters")
inst = InstrumentInfo(imo, "$rel/<telescope>/instrument_info")
ch = FreqChannelInfo(imo, "$rel/<telescope>/<channel>/channel_info")
det = DetectorInfo(imo, "$rel/<telescope>/<channel>/<detector>/detector_info")
```

[`query`](@ref) accepts the same identifiers as `Imo.query` in litebird_sim: release paths (the `/releases`
prefix is optional), UUID strings, and `/data_files/`, `/quantities/` or `/entities/` followed by a UUID. If
a release holds several versions of a quantity, the most recently uploaded one is used. Only uncompressed
`schema.json` files are read; decompress `schema.json.gz` first.

The constructors mirror litebird_sim's `from_imo` / `from_dict`: angles stored in degrees are converted to
radians, `FreqChannelInfo` derives missing NETs and detector names, and unknown keys are ignored. A
time-dependent detector quaternion is a dictionary with `start_time` (seconds), `sampling_rate_hz` and
`quats`; calendar start times are not supported. Unlike litebird_sim, `pixel` and `squid` are stored as
integers and no default bandpass is synthesised.

## Lists of detectors

[`detector_list_from_parameters`](@ref) builds detectors from dictionaries with litebird_sim's keys:

```julia
dets = detector_list_from_parameters(imo, [
    Dict("channel_info_obj" => "$rel/<telescope>/<channel>/channel_info",
         "num_of_detectors_from_channel" => 4),
    Dict("detector_info_obj" => "$rel/<telescope>/<channel>/<detector>/detector_info",
         "fwhm_arcmin" => 30.0),                        # fields override the IMO values
    Dict("name" => "mock", "sampling_rate_hz" => 19.0),
])
```

## Parameter files

[`run_h_maps_pipeline`](@ref) reads a TOML file whose sections follow litebird_sim's names, computes the
h-maps of the detectors and saves them with [`save_h_maps`](@ref):

```toml
[simulation]
base_path = "output"
start_time = 0.0
duration_s = "1 year"            # seconds, or "N years/days/hours/min/s"

[imo]
flatfile_location = "/path/to/imo"

[scanning_strategy]
imo_url = "/releases/<tag>/satellite/scanning_parameters"
# delta_time_s = 2.0             # optional: litebird_sim's interpolation

[instrument]
imo_url = "/releases/<tag>/satellite/<telescope>/instrument_info"

[hwp]
ang_speed_radpsec = 4.817        # 2π · rpm / 60

[[detectors]]
channel_info_obj = "/releases/<tag>/satellite/<telescope>/<channel>/channel_info"

[h_maps]
nside = 128
n_m_couples = [[0, 0], [2, 0], [4, 0], [2, -4], [2, 4]]
output_coordinate_system = "Galactic"
output_directory = "h_maps"      # relative to base_path
```

```julia
res = run_h_maps_pipeline("params.toml")
```

The scanning strategy and the instrument can also be given explicitly (`spin_sun_angle_rad`,
`precession_rate_hz`, `spin_rate_hz`; `spin_boresight_angle_rad`, `spin_rotangle_rad`,
`boresight_rotangle_rad`) instead of `imo_url`.
