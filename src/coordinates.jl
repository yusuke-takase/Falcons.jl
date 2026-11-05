"""
    CoordinateSystem

Coordinate system of pointings and maps: `Ecliptic` or `Galactic`, as
`CoordinateSystem` in litebird_sim.
"""
@enum CoordinateSystem Ecliptic Galactic

@doc "Ecliptic coordinates (`BarycentricMeanEcliptic` in astropy), the default of `get_pointings`." Ecliptic
@doc "Galactic coordinates, the default of `get_hitmap`, `get_pointing_pixels` and `make_h_maps`." Galactic

"""
Rotation matrix from Ecliptic to Galactic coordinates, `v_gal =
ECL_TO_GAL_ROT_MATRIX * v_ecl`. Same values as litebird_sim, which computes the
matrix with astropy (`BarycentricMeanEcliptic` → `galactic`, astropy 8.0.2).
"""
const ECL_TO_GAL_ROT_MATRIX = @SMatrix [
	-0.05487563486628303 -0.993821352389173 -0.09647685853808491
	0.49410952885789666 -0.1109909720383276 0.8622857865103742
	-0.8676660868040481 -0.0003516551470140186 0.497147300252961
]

"""
Quaternion of the rotation from Ecliptic to Galactic coordinates, as
`ECL_TO_GAL_QUAT` in litebird_sim. Multiplied on the left of the
boresight-to-Ecliptic quaternions, it gives pointings in Galactic coordinates.
"""
const ECL_TO_GAL_QUAT = Quaternion(-0.37382082503489683, 0.33419207148546864, 0.6447894825362716, 0.5769056884249691)

"The Ecliptic north pole in Galactic coordinates."
const NORTH_POLE_VEC = ECL_TO_GAL_ROT_MATRIX[:, 3]

"Quaternion that rotates Ecliptic coordinates into `coord`."
coordinate_quat(coord::CoordinateSystem) = coord == Galactic ? ECL_TO_GAL_QUAT : one(Quaternion{Float64})

"""
    rotate_vector_e2g!(v) -> v

Rotate the 3-vector `v` in place from Ecliptic to Galactic coordinates.
"""
function rotate_vector_e2g!(v::AbstractVector{<:AbstractFloat})
	v .= ECL_TO_GAL_ROT_MATRIX * SVector{3}(v[1], v[2], v[3])
	return v
end

function _ang2galvec_one_sample(theta, phi)
	st = sin(theta)
	return ECL_TO_GAL_ROT_MATRIX * SVector(st * cos(phi), st * sin(phi), cos(theta))
end

function _rotate_coordinates_and_pol_e2g_for_one_sample(theta_ecl, phi_ecl, pol_angle_ecl)
	vec = _ang2galvec_one_sample(theta_ecl, phi_ecl)
	theta_gal = atan(sqrt(vec[1]^2 + vec[2]^2), vec[3])
	phi_gal = _wrap_phi(atan(vec[2], vec[1]))
	sinalpha = NORTH_POLE_VEC[1] * vec[2] - NORTH_POLE_VEC[2] * vec[1]
	cosalpha = NORTH_POLE_VEC[3] - vec[3] * (NORTH_POLE_VEC ⋅ vec)
	pol_angle_gal = pol_angle_ecl + atan(sinalpha, cosalpha)
	return (theta_gal, phi_gal, pol_angle_gal)
end

"""
    rotate_coordinates_e2g!(θ, φ, ψ)

Rotate arrays of pointings `(θ, φ)` and orientation angles `ψ` (rad) in place
from Ecliptic to Galactic coordinates, as `rotate_coordinates_e2g` in
litebird_sim (which follows healpy). `ψ` is not wrapped afterwards.
"""
function rotate_coordinates_e2g!(θ, φ, ψ)
	@inbounds for i in eachindex(θ, φ, ψ)
		θ[i], φ[i], ψ[i] = _rotate_coordinates_and_pol_e2g_for_one_sample(θ[i], φ[i], ψ[i])
	end
	return (θ, φ, ψ)
end
