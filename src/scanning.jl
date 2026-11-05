"""
Angular velocity of the anti-Sun direction (rad/s): one revolution per Julian
year of 365.25 days, as `YEARLY_OMEGA_SPIN_HZ` in litebird_sim.
"""
const YEARLY_OMEGA_SPIN_HZ = 2π / (365.25 * 86400.0)

"""
    SpinningScanningStrategy(; spin_sun_angle_rad, precession_rate_hz, spin_rate_hz)

Scanning strategy of a spinning spacecraft whose spin axis precesses around the
anti-Sun direction, as `SpinningScanningStrategy` in litebird_sim.

- `spin_sun_angle_rad`: angle α between the spin axis and the anti-Sun direction.
- `precession_rate_hz`: precession frequency (rotations per second).
- `spin_rate_hz`: spin frequency (rotations per second).

The spin-to-Ecliptic rotation at the absolute time `t` (s) is

    Rz(Ω t) · Rx(2π f_prec t) · Ry(π/2 − α) · Rz(2π f_spin t)

with Ω = [`YEARLY_OMEGA_SPIN_HZ`](@ref); see [`spin_to_ecliptic`](@ref).
"""
Base.@kwdef struct SpinningScanningStrategy
	spin_sun_angle_rad::Float64
	precession_rate_hz::Float64
	spin_rate_hz::Float64
end

"""
    SpinningScanningStrategy(imo::Imo, url)

Scanning strategy from an IMO data file with `spin_sun_angle_deg`,
`precession_period_min` and `spin_rate_rpm`, as
`SpinningScanningStrategy.from_imo` in litebird_sim.
"""
function SpinningScanningStrategy(imo::Imo, url::AbstractString)
	m = _metadata(imo, url)
	return SpinningScanningStrategy(
		spin_sun_angle_rad = deg2rad(Float64(m["spin_sun_angle_deg"])),
		precession_rate_hz = 1.0 / (60.0 * m["precession_period_min"]),
		spin_rate_hz = m["spin_rate_rpm"] / 60.0,
	)
end

function Base.show(io::IO, ::MIME"text/plain", s::SpinningScanningStrategy)
	println(io, "SpinningScanningStrategy")
	println(io, "  spin_sun_angle_rad : ", s.spin_sun_angle_rad, " (", rad2deg(s.spin_sun_angle_rad), "°)")
	println(io, "  precession_rate_hz : ", s.precession_rate_hz, " (period ", 1 / (60 * s.precession_rate_hz), " min)")
	print(io, "  spin_rate_hz       : ", s.spin_rate_hz, " (", 60 * s.spin_rate_hz, " rpm)")
end

"""
    spin_to_ecliptic(scan::SpinningScanningStrategy, t) -> Quaternion

Rotation from the spin-axis frame to Ecliptic coordinates at the absolute time
`t` (s). Same operations, in the same order, as `spin_to_ecliptic` in
litebird_sim with a floating-point start time.
"""
@inline function spin_to_ecliptic(scan::SpinningScanningStrategy, t::Real)
	q = quat_rotation_z(2π * scan.spin_rate_hz * t)
	q = quat_rotation_y(π / 2 - scan.spin_sun_angle_rad) * q
	q = quat_rotation_x(2π * scan.precession_rate_hz * t) * q
	q = quat_rotation_z(YEARLY_OMEGA_SPIN_HZ * t) * q
	return q
end

@inline _clip_sincos(x) = min(max(x, -1), 1)

"""
    orientation_angle(θ, φ, ordir) -> ψ

Orientation angle (rad, in (−π, π]) of the direction `ordir` at the sky
position `(θ, φ)`, measured from the local East (`e_φ`) towards the local North
(`−e_θ`), as `orientation_angle` in litebird_sim.
"""
@inline function orientation_angle(θ, φ, ordir)
	cos_psi = _clip_sincos(-sin(φ) * ordir[1] + cos(φ) * ordir[2])
	sin_psi = _clip_sincos((-cos(θ) * cos(φ) * ordir[1]) + (-cos(θ) * sin(φ) * ordir[2]) + (sin(θ) * ordir[3]))
	return atan(sin_psi, cos_psi)
end

@inline function _wrap_phi(φ)
	φ < 0 && (φ += 2π)
	φ >= 2π && (φ = zero(φ))
	return φ
end

"""
    compute_pointing_and_orientation(q::Quaternion) -> (θ, φ, ψ)

Sky direction `(θ, φ)` of the z axis rotated by `q`, with `φ ∈ [0, 2π)`, and the
orientation angle `ψ` of the rotated x axis, as
`compute_pointing_and_orientation` in litebird_sim.
"""
@inline function compute_pointing_and_orientation(q::Quaternion)
	v = rotate_z_vector(q)
	θ = atan(sqrt(v[1]^2 + v[2]^2), v[3])
	φ = atan(v[2], v[1])
	ψ = orientation_angle(θ, φ, rotate_x_vector(q))
	return (θ, _wrap_phi(φ), ψ)
end

# name => (α deg, β deg, precession period min, spin period min)
const SATELLITE_PRESETS = (
	WMAP = (70.0, 22.5, 60.0, 129 / 60),
	Planck = (7.5, 85.0, 24 * 60 * (365.25 / 2), 1.0),
	PICO = (26.0, 69.0, 10 * 60.0, 1.0),
	CORE = (30.0, 65.0, 4 * 24 * 60.0, 2.0),
	EPIC = (45.0, 55.0, 3.2 * 60, 1.0),
	LiteBIRD = (45.0, 50.0, 192.348, 20.0),
)

"""
    satellite_presets() -> Vector{Symbol}

Names accepted by [`satellite_preset`](@ref).
"""
satellite_presets() = collect(keys(SATELLITE_PRESETS))

"""
    satellite_preset(name::Symbol) -> (SpinningScanningStrategy, InstrumentInfo)

Scanning strategy and instrument geometry (spin-boresight angle β) of a known
satellite design: `:WMAP`, `:Planck`, `:PICO`, `:CORE`, `:EPIC` or `:LiteBIRD`.
"""
function satellite_preset(name::Symbol)
	haskey(SATELLITE_PRESETS, name) ||
		throw(ArgumentError("unknown satellite :$name; choose one of $(join(":" .* string.(satellite_presets()), ", "))"))
	α, β, prec_period_min, spin_period_min = SATELLITE_PRESETS[name]
	scan = SpinningScanningStrategy(
		spin_sun_angle_rad = deg2rad(α),
		precession_rate_hz = 1 / (60 * prec_period_min),
		spin_rate_hz = 1 / (60 * spin_period_min),
	)
	return scan, InstrumentInfo(name = string(name), spin_boresight_angle_rad = deg2rad(β))
end
