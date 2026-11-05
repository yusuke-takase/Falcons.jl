# Time-dependent rotations and external attitude

Any rotation of the chain `q(t) = A(t) · S(t) · bore2spin · d(t)` (see [Conventions](@ref)) except
`bore2spin` can be a time series. Time series are [`RotQuaternion`](@ref)s: quaternions sampled at the
absolute times `start_time_s + k / sampling_rate_hz` and interpolated with [`slerp`](@ref).

```julia
using Falcons

q0 = detector_quat(deg2rad(2), deg2rad(120))
t = 0.0:1.0:4000.0
wobble = [quat_rotation_x(deg2rad(1 / 60) * sin(2π * tk / 97)) * q0 for tk in t]
d = RotQuaternion(wobble; start_time_s = 0.0, sampling_rate_hz = 1.0)
```

A series must cover every time at which it is evaluated: otherwise the computation stops with an
`ArgumentError` (litebird_sim silently wraps around instead).

## Detector quaternions

Pass a time series as `DetectorInfo(quat = d)`. In the IMO the same information is a dictionary with
`start_time` (seconds), `sampling_rate_hz` and `quats` (see [Instrument database and parameter files](@ref)).

## External attitude

The first argument of an [`Observation`](@ref) is the attitude `A(t)`, the rotation from the spin-axis frame
to Ecliptic coordinates. Instead of a [`SpinningScanningStrategy`](@ref) it can be a time series, for
instance a measured or simulated attitude:

```julia
# quats::Vector{Quaternion{Float64}}: spin-axis frame → Ecliptic, sampled every second from t0
attitude = RotQuaternion(quats; start_time_s = t0, sampling_rate_hz = 1.0)
obs = Observation(attitude, inst, dets; start_time_s = t0, duration_s = 3600.0)
```

Falcons does not read attitude files: convert the data to this frame before building the
`RotQuaternion`. Constant frame changes can be applied to the whole series with a product, e.g.
`inertial_to_ecl * RotQuaternion(quats; ...)` for a rotation on the left (from the data's inertial frame to
Ecliptic) or `RotQuaternion(quats; ...) * spin_to_body` on the right (from the spin-axis frame to the data's
body frame). The angles of the [`InstrumentInfo`](@ref) are then defined with respect to that spin-axis frame.

## Perturbations of the nominal attitude

To keep the analytic scan and add a small rotation in the spin-axis frame (for example a pointing error or
a wobble), pass it as `spin_perturbation`:

```julia
obs = Observation(scan, inst, dets; start_time_s = 0.0, duration_s = 3600.0,
                  spin_perturbation = RotQuaternion(dq; start_time_s = 0.0, sampling_rate_hz = 10.0))
```

This is the position where litebird_sim's `SpacecraftCoord` adds its disturbances.

## Exact and interpolated pointings

By default every sample is computed exactly: the analytic attitude is evaluated at the sample time and the
time series are interpolated to it. With `Observation(...; delta_time_s = Δ)` Falcons reproduces
litebird_sim instead: the boresight rotation `A · S · bore2spin` is computed on anchors every `Δ` seconds
(counted from `start_time_simulation_s`) and interpolated, and time-dependent detector quaternions are
interpolated separately. The anchors are computed for one block of samples at a time, so their memory does
not grow with the length of the observation. Both modes are available on the CPU and on the GPU.
