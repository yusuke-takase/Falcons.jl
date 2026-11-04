# Conventions

Falcons uses the conventions of [litebird_sim](https://github.com/litebird/litebird_sim). The test suite
compares Falcons with litebird_sim on reference configurations (`test/reference/`), so the statements below
are checked numerically.

## Quaternions

[`Quaternion`](@ref)s are stored as `(x, y, z, w)`, products follow the Hamilton convention and rotations
are active: a vector `v` is rotated to `q v q*` ([`rotate_vector`](@ref)). `q1 * q2` applies `q2` first.
[`quat_rotation`](@ref)`(θ, axis)` and `quat_rotation_x/y/z(θ)` build right-handed rotations by `θ` radians.

## Frames and the rotation chain

- **Detector frame**: the beam points along `+z`, the polarisation reference along `+x`.
- **Boresight frame**: a detector is placed with its quaternion `d` (detector → boresight), see
  [`detector_quat`](@ref).
- **Spin-axis frame**: the boresight frame is rotated by
  `bore2spin = Rz(spin_rotangle) · Ry(spin_boresight_angle) · Rz(boresight_rotangle)`
  ([`InstrumentInfo`](@ref)).
- **Ecliptic frame**: the spin-axis frame is rotated by the attitude
  `A(t) = Rz(Ω t) · Rx(2π f_prec t) · Ry(π/2 − α) · Rz(2π f_spin t)`
  ([`SpinningScanningStrategy`](@ref), [`spin_to_ecliptic`](@ref)), with α the spin–anti-Sun angle and
  `Ω = 2π / (365.25 days)`. At `t = 0` the anti-Sun direction is `+x` and the Ecliptic north pole is `+z`.

The quaternion of a detector at the time `t` is therefore

```
q(t) = A(t) · S(t) · bore2spin · d(t)
```

where `S(t)` is an optional perturbation in the spin-axis frame and `A`, `S` and `d` can also be time
series (see [Time-dependent rotations and external attitude](@ref)). For Galactic outputs the Ecliptic-to-
Galactic rotation is multiplied on the left; its values are litebird_sim's (astropy,
`BarycentricMeanEcliptic` → `galactic`).

**Times are absolute.** An [`Observation`](@ref) samples `t_k = start_time_s + k / sampling_rate_hz`, and
the spin, precession and anti-Sun phases are evaluated at `t_k` itself, as litebird_sim does for a
floating-point start time. Calendar dates and ephemerides are not supported.

## Pointing angles

For every sample [`get_pointings`](@ref) returns `(θ, φ, ψ)`:

- `θ ∈ [0, π]`: colatitude of the beam direction `q ẑ q*`;
- `φ ∈ [0, 2π)`: its longitude;
- `ψ ∈ (−π, π]`: orientation of `q x̂ q*`, measured from the local East (`e_φ`) towards the local North
  (`−e_θ`).

The detector polarisation angle (`DetectorInfo.pol_angle_rad`) is **not** included in ψ.

## Half-wave plate

An [`IdealHWP`](@ref) has the mechanical angle

```
χ(t) = start_angle_rad + ang_speed_radpsec · (t − start_time_simulation_s)
```

wrapped to [−π, π]. It is measured from the start of the simulation, so observations split with
[`split_observation`](@ref) keep a continuous phase. The angle is evaluated in closed form, while
litebird_sim accumulates it with a recurrence; the two differ by rounding only (about 1e-12 rad after ten
minutes at 19 Hz).

In litebird_sim's signal model the polarisation angle seen by a detector behind an ideal HWP is
`2χ − pol_angle_rad + ψ`.

## h-maps

[`make_h_maps`](@ref) computes, in every pixel `p` and for integer spin orders `(n, m)`,

```
h_{n,m}(p) = (1 / N_p) Σ_{samples in p} exp(+i (n ψ + m χ))
```

with `N_p` the number of hits, as `make_h_maps` in litebird_sim. For `(n, m) = (0, 0)` the real part holds
`N_p` instead, and pixels without hits are 0. Observations of the same detector (matched by name) are
accumulated before normalising. Files written by [`save_h_maps`](@ref) use litebird_sim's layout
(`h_maps_det_<name>.h5`, one group `"n,m"` with `Re` and `Im` per couple).

## Pixelisation

Maps use HEALPix with 1-based pixel indices (Healpix.jl's convention; subtract 1 for healpy). Hit maps,
pixel indices and h-maps default to Galactic coordinates and RING ordering, as in litebird_sim; pointings
default to Ecliptic coordinates.
