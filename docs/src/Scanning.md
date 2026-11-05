## Scanning strategy

A scan is described with the types of [litebird_sim](https://github.com/litebird/litebird_sim), with the
same fields and units (angles in radians, rates in Hz):

```julia
using Falcons

scan = SpinningScanningStrategy(
    spin_sun_angle_rad = deg2rad(45),
    precession_rate_hz = 1 / (192.348 * 60),
    spin_rate_hz = 0.05 / 60,
)
inst = InstrumentInfo(spin_boresight_angle_rad = deg2rad(50))
```

`satellite_preset` returns the scanning strategy and the instrument geometry of a few known designs:

```julia
scan, inst = satellite_preset(:LiteBIRD)
```

Detectors are `DetectorInfo`s. `detector_quat(θ, φ; ψ)` places a detector at the focal-plane position
`(θ, φ)` (rad, measured from the boresight) with its polarisation reference turned by `ψ`:

```julia
dets = [
    DetectorInfo(name = "boresight", sampling_rate_hz = 19.0),
    DetectorInfo(name = "offset", sampling_rate_hz = 19.0, quat = detector_quat(deg2rad(5), deg2rad(30))),
]
```

## Generate pointings

An `Observation` groups the scan, the instrument and the detectors over a time interval. Times are absolute
seconds; an `IdealHWP` adds a rotating half-wave plate.

```julia
obs = Observation(scan, inst, dets; start_time_s = 0.0, duration_s = 3600.0,
                  hwp = IdealHWP(ang_speed_radpsec = 46 * 2π / 60))

pointings, hwp_angle = get_pointings(obs)    # pointings[:, k, i] = (θ, φ, ψ) of detector i at sample k
px = get_pointing_pixels(obs; nside = 128)  # px.pixels, px.psi and px.hwp_angle, without storing (θ, φ)
```

`θ` and `φ` are the colatitude and longitude of the beam (Ecliptic coordinates by default for
`get_pointings`, Galactic for the pixel functions) and `ψ` is the orientation angle, defined as in
litebird_sim.

## Hit maps and h-maps

```julia
hits = get_hitmap(obs; nside = 128)
res = make_h_maps(obs; nside = 128, n_m_couples = nm_grid(0:4, (-4, 0, 4)))
h20 = h_map(res, "boresight", 2, 0)          # h_{2,0} of one detector
save_h_maps(res, "h_maps")                    # one HDF5 file per detector, litebird_sim's layout
```

The h-maps are `h_{n,m}(p) = (1 / N_p) Σ exp(i (n ψ + m χ))` over the samples in pixel `p`, with `χ` the HWP
angle; the `(0, 0)` couple holds the hit counts.
