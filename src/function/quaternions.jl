"""
    Quaternion{T}(x, y, z, w)

Quaternion with vector part `(x, y, z)` and scalar part `w`, stored in the same
`(x, y, z, w)` order as litebird_sim and as the `quat` entries of the IMO.

Products follow the Hamilton convention and rotations are active: the vector
`v` is rotated to `q v q*` (see [`rotate`](@ref)). The type is `isbits`, so it can
be passed to GPU kernels.
"""
struct Quaternion{T <: Real}
	x::T
	y::T
	z::T
	w::T
end

Quaternion(x::Real, y::Real, z::Real, w::Real) = Quaternion(promote(float(x), y, z, w)...)

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

"""
    vector_part(q::Quaternion) -> SVector{3}
"""
vector_part(q::Quaternion) = SVector(q.x, q.y, q.z)

"""
    rotation_quat(axis, angle) -> Quaternion

Quaternion of the right-handed rotation by `angle` (rad) around the unit vector
`axis`.
"""
@inline function rotation_quat(axis, angle)
	s, c = sincos(angle / 2)
	return Quaternion(axis[1] * s, axis[2] * s, axis[3] * s, c)
end

"""
    rotate(q::Quaternion, v) -> SVector{3}

Rotate the 3-vector `v` by the unit quaternion `q`, i.e. return the vector part
of `q v q*`.
"""
@inline function rotate(q::Quaternion, v)
	# v' = v + 2 w (u × v) + 2 u × (u × v), with u the vector part of q.
	u = SVector(q.x, q.y, q.z)
	vv = SVector{3}(v[1], v[2], v[3])
	t = 2 * cross(u, vv)
	return vv + q.w * t + cross(u, t)
end
