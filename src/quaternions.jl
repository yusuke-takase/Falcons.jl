"""
    Quaternion{T}(x, y, z, w)

Quaternion with vector part `(x, y, z)` and scalar part `w`, stored in the same
`(x, y, z, w)` order as litebird_sim and as the `quat` entries of the IMO.

Products follow the Hamilton convention and rotations are active: the vector
`v` is rotated to `q v q*` (see [`rotate_vector`](@ref)). The type is `isbits`,
so it can be passed to GPU kernels.

`Quaternion(x, y, z, w)` promotes the components to a common floating-point
type, so `Quaternion(0, 0, 0, 1)` is a `Quaternion{Float64}`.
"""
struct Quaternion{T <: Real}
	x::T
	y::T
	z::T
	w::T

	Quaternion{T}(x, y, z, w) where {T <: Real} = new{T}(x, y, z, w)
end

function Quaternion(x::Real, y::Real, z::Real, w::Real)
	c = promote(float(x), float(y), float(z), float(w))
	return Quaternion{eltype(c)}(c...)
end

"""
    Quaternion(v::AbstractVector)

Build a quaternion from a 4-element vector in `(x, y, z, w)` order.
"""
function Quaternion(v::AbstractVector)
	length(v) == 4 || throw(ArgumentError("a quaternion needs 4 components (x, y, z, w), got $(length(v))"))
	return Quaternion(v[1], v[2], v[3], v[4])
end

Base.one(::Type{Quaternion{T}}) where {T} = Quaternion{T}(zero(T), zero(T), zero(T), one(T))
Base.one(::Quaternion{T}) where {T} = one(Quaternion{T})

function Base.:*(a::Quaternion, b::Quaternion)
	return Quaternion(
		a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
		a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
		a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
		a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
	)
end

Base.:-(q::Quaternion) = Quaternion(-q.x, -q.y, -q.z, -q.w)
Base.float(q::Quaternion) = Quaternion(float(q.x), float(q.y), float(q.z), float(q.w))
Base.conj(q::Quaternion) = Quaternion(-q.x, -q.y, -q.z, q.w)
Base.abs2(q::Quaternion) = q.x^2 + q.y^2 + q.z^2 + q.w^2
Base.abs(q::Quaternion) = sqrt(abs2(q))

function Base.inv(q::Quaternion)
	n = abs2(q)
	return Quaternion(-q.x / n, -q.y / n, -q.z / n, q.w / n)
end

Base.isapprox(a::Quaternion, b::Quaternion; kw...) = isapprox(components(a), components(b); kw...)

"""
    components(q::Quaternion) -> SVector{4}

The components of `q` in `(x, y, z, w)` order.
"""
components(q::Quaternion) = SVector(q.x, q.y, q.z, q.w)

"Dot product of the components of two quaternions."
quat_dot(a::Quaternion, b::Quaternion) = a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w

"""
    normalize(q::Quaternion) -> Quaternion

`q` divided by its norm.
"""
function LinearAlgebra.normalize(q::Quaternion)
	n = abs(q)
	return Quaternion(q.x / n, q.y / n, q.z / n, q.w / n)
end

"""
    quat_rotation(θ, axis) -> Quaternion

Quaternion of the right-handed rotation by `θ` (rad) around the unit vector
`axis`. Same as `quat_rotation` in litebird_sim.
"""
@inline function quat_rotation(θ, axis)
	s, c = sincos(θ / 2)
	return Quaternion(axis[1] * s, axis[2] * s, axis[3] * s, c)
end

"Quaternion of the rotation by `θ` (rad) around the x axis."
@inline quat_rotation_x(θ) = (s = sin(θ / 2); Quaternion(s, zero(s), zero(s), cos(θ / 2)))
"Quaternion of the rotation by `θ` (rad) around the y axis."
@inline quat_rotation_y(θ) = (s = sin(θ / 2); Quaternion(zero(s), s, zero(s), cos(θ / 2)))
"Quaternion of the rotation by `θ` (rad) around the z axis."
@inline quat_rotation_z(θ) = (s = sin(θ / 2); Quaternion(zero(s), zero(s), s, cos(θ / 2)))

"""
    rotate_vector(q::Quaternion, v) -> SVector{3}

Rotate the 3-vector `v` by the unit quaternion `q`, i.e. return the vector part
of `q v q*`. Same as `rotate_vector` in litebird_sim.
"""
@inline function rotate_vector(q::Quaternion, v)
	# v' = v + 2 w (u × v) + 2 u × (u × v), with u the vector part of q.
	u = SVector(q.x, q.y, q.z)
	vv = SVector{3}(v[1], v[2], v[3])
	t = 2 * cross(u, vv)
	return vv + q.w * t + cross(u, t)
end

"The x axis rotated by the unit quaternion `q` (litebird_sim's `rotate_x_vector`)."
@inline function rotate_x_vector(q::Quaternion)
	return SVector(1 - 2 * (q.y * q.y + q.z * q.z), 2 * (q.x * q.y + q.w * q.z), 2 * (q.x * q.z - q.w * q.y))
end

"The z axis rotated by the unit quaternion `q` (litebird_sim's `rotate_z_vector`)."
@inline function rotate_z_vector(q::Quaternion)
	return SVector(2 * (q.w * q.y + q.x * q.z), 2 * (q.y * q.z - q.w * q.x), 1 - 2 * (q.x * q.x + q.y * q.y))
end

"""
    slerp(a::Quaternion, b::Quaternion, t) -> Quaternion

Spherical linear interpolation between the unit quaternions `a` (`t = 0`) and
`b` (`t = 1`) along the shortest path.
"""
@inline function slerp(a::Quaternion, b::Quaternion, t::Real)
	d = quat_dot(a, b)
	if d < 0
		b = -b
		d = -d
	end
	if d > 1 - 1e-12
		# Nearly identical rotations: linear interpolation is exact to rounding.
		return normalize(Quaternion(
			a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), a.z + t * (b.z - a.z), a.w + t * (b.w - a.w)))
	end
	ω = acos(d)
	s = sin(ω)
	ca, cb = sin((1 - t) * ω) / s, sin(t * ω) / s
	return Quaternion(ca * a.x + cb * b.x, ca * a.y + cb * b.y, ca * a.z + cb * b.z, ca * a.w + cb * b.w)
end

"""
    RotQuaternion(q::Quaternion)
    RotQuaternion(quats::AbstractVector{<:Quaternion}; start_time_s, sampling_rate_hz)
    RotQuaternion(quats::AbstractMatrix; start_time_s, sampling_rate_hz)

A rotation that is either constant or sampled in time, as `RotQuaternion` in
litebird_sim.

A time series holds the quaternions `quats[k]` at the absolute times
`start_time_s + (k - 1) / sampling_rate_hz` (s); values in between are obtained
by spherical linear interpolation ([`slerp`](@ref)). Asking for a time outside
the covered range throws an `ArgumentError`. A matrix argument must have 4
rows, one quaternion `(x, y, z, w)` per column. Quaternions are normalised on
construction.

`quat * r` and `r * quat` with a constant `Quaternion` apply it to every sample;
the product of two time series requires the same time grid.
"""
struct RotQuaternion{T <: Real}
	quats::Vector{Quaternion{T}}
	start_time_s::Union{Nothing, Float64}
	sampling_rate_hz::Union{Nothing, Float64}

	function RotQuaternion{T}(quats::Vector{Quaternion{T}}, start_time_s, sampling_rate_hz) where {T}
		isempty(quats) && throw(ArgumentError("a RotQuaternion needs at least one quaternion"))
		if start_time_s === nothing
			length(quats) == 1 || throw(ArgumentError("a constant RotQuaternion holds exactly one quaternion; pass start_time_s and sampling_rate_hz for a time series"))
			sampling_rate_hz === nothing || throw(ArgumentError("a constant RotQuaternion has no sampling rate"))
		else
			length(quats) >= 2 || throw(ArgumentError("a time-dependent RotQuaternion needs at least two quaternions"))
			(sampling_rate_hz !== nothing && sampling_rate_hz > 0) || throw(ArgumentError("sampling_rate_hz must be positive"))
		end
		return new{T}(normalize.(quats), start_time_s === nothing ? nothing : Float64(start_time_s),
			sampling_rate_hz === nothing ? nothing : Float64(sampling_rate_hz))
	end
end

# The quaternions are stored with a floating-point element type, so that they
# can hold their normalised values.
_rotquat(quats::Vector{Quaternion{T}}, start_time_s, sampling_rate_hz) where {T} =
	RotQuaternion{T}(quats, start_time_s, sampling_rate_hz)

RotQuaternion(q::Quaternion) = _rotquat([float(q)], nothing, nothing)
RotQuaternion(q::RotQuaternion) = q
RotQuaternion(v::AbstractVector{<:Real}) = RotQuaternion(Quaternion(v))

function RotQuaternion(quats::AbstractVector{<:Quaternion}; start_time_s::Real, sampling_rate_hz::Real)
	return _rotquat([float(q) for q in quats], start_time_s, sampling_rate_hz)
end

function RotQuaternion(quats::AbstractMatrix{<:Real}; start_time_s::Real, sampling_rate_hz::Real)
	size(quats, 1) == 4 || throw(ArgumentError("the matrix must have 4 rows (x, y, z, w), got size $(size(quats))"))
	return RotQuaternion([Quaternion(view(quats, :, k)) for k in axes(quats, 2)]; start_time_s, sampling_rate_hz)
end

"`true` if `r` is constant in time."
is_constant(r::RotQuaternion) = r.start_time_s === nothing

Base.length(r::RotQuaternion) = length(r.quats)

"""
    time_range(r::RotQuaternion) -> (first_s, last_s)

First and last sample time of a time-dependent rotation.
"""
function time_range(r::RotQuaternion)
	is_constant(r) && return (-Inf, Inf)
	return (r.start_time_s, r.start_time_s + (length(r.quats) - 1) / r.sampling_rate_hz)
end

"""
    check_covers(r::RotQuaternion, t_first, t_last, what)

Throw an `ArgumentError` unless the time series `r` covers `[t_first, t_last]`.
"""
function check_covers(r::RotQuaternion, t_first, t_last, what::AbstractString)
	is_constant(r) && return nothing
	r_first, r_last = time_range(r)
	tol = 1e-9 / r.sampling_rate_hz
	if t_first < r_first - tol || t_last > r_last + tol
		throw(ArgumentError("$what covers $(r_first) – $(r_last) s, but $(t_first) – $(t_last) s are needed"))
	end
	return nothing
end

"""
    slerp(r::RotQuaternion, t) -> Quaternion

The rotation `r` at the absolute time `t` (s). Throws an `ArgumentError` when
`t` is outside the range covered by a time series.
"""
function slerp(r::RotQuaternion, t::Real)
	is_constant(r) && return r.quats[1]
	fi = (t - r.start_time_s) * r.sampling_rate_hz
	n = length(r.quats)
	if !(-1e-9 <= fi <= n - 1 + 1e-9)
		first_s, last_s = time_range(r)
		throw(ArgumentError("time $t s is outside the range $first_s – $last_s s of the RotQuaternion"))
	end
	i = clamp(floor(Int, fi), 0, n - 2)
	return slerp(r.quats[i+1], r.quats[i+2], fi - i)
end

Base.:*(q::Quaternion, r::RotQuaternion) = _rotquat_like(r, [q * p for p in r.quats])
Base.:*(r::RotQuaternion, q::Quaternion) = _rotquat_like(r, [p * q for p in r.quats])

function Base.:*(a::RotQuaternion, b::RotQuaternion)
	is_constant(a) && return a.quats[1] * b
	is_constant(b) && return a * b.quats[1]
	if length(a) != length(b) || !isapprox(a.start_time_s, b.start_time_s; atol = 1e-9 / a.sampling_rate_hz, rtol = 0) ||
	   !isapprox(a.sampling_rate_hz, b.sampling_rate_hz; rtol = 1e-12)
		throw(ArgumentError("the two time series have different time grids; resample one with `slerp` first"))
	end
	return _rotquat_like(a, [p * q for (p, q) in zip(a.quats, b.quats)])
end

function _rotquat_like(r::RotQuaternion, quats::Vector{<:Quaternion})
	is_constant(r) && return RotQuaternion(quats[1])
	return RotQuaternion(quats; start_time_s = r.start_time_s, sampling_rate_hz = r.sampling_rate_hz)
end
