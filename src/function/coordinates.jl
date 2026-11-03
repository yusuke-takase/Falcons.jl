"""
Rotation matrix from Ecliptic (J2000) to Galactic coordinates: `v_gal =
ECL_TO_GAL_ROT_MATRIX * v_ecl`. The columns are the Ecliptic axes expressed in
Galactic coordinates, as computed with SkyCoords.jl (`EclipticCoords{2000}` →
`GalCoords`).
"""
const ECL_TO_GAL_ROT_MATRIX = @SMatrix [
	-0.05487553939571462 -0.9938213632550472 -0.09647680091042089
	0.49410945362336356 -0.11099087460650042 0.8622858421626236
	-0.8676661356858699 -0.0003516986832161391 0.4971472149092117
]

"The Ecliptic north pole in Galactic coordinates."
const NORTH_POLE_VEC = ECL_TO_GAL_ROT_MATRIX[:, 3]

"""
    rotate_vec_ecl_to_gal_coordinates!(v) -> v

Rotate the 3-vector `v` in place from Ecliptic to Galactic coordinates.
"""
function rotate_vec_ecl_to_gal_coordinates!(v::AbstractVector{<:AbstractFloat})
	v .= ECL_TO_GAL_ROT_MATRIX * SVector{3}(v[1], v[2], v[3])
	return v
end

"""
    _ang2galvec_one_sample(theta, phi) -> SVector{3}

Unit vector in Galactic coordinates of the Ecliptic direction `(theta, phi)`
(colatitude and longitude, rad). Follows `healpy.rotator`.
"""
function _ang2galvec_one_sample(theta, phi)
	st = sin(theta)
	return ECL_TO_GAL_ROT_MATRIX * SVector(st * cos(phi), st * sin(phi), cos(theta))
end

"""
    _rotate_coordinates_and_pol_e2g_for_one_sample(theta, phi, psi) -> (θ, φ, ψ)

Rotate a pointing `(theta, phi)` and orientation angle `psi` (rad) from
Ecliptic to Galactic coordinates. Follows litebird_sim and healpy.
"""
function _rotate_coordinates_and_pol_e2g_for_one_sample(theta_ecl, phi_ecl, pol_angle_ecl)
	vec = _ang2galvec_one_sample(theta_ecl, phi_ecl)
	theta_gal, phi_gal = vec2ang(vec[1], vec[2], vec[3])
	sinalpha = NORTH_POLE_VEC[1] * vec[2] - NORTH_POLE_VEC[2] * vec[1]
	cosalpha = NORTH_POLE_VEC[3] - vec[3] * (NORTH_POLE_VEC ⋅ vec)
	pol_angle_gal = pol_angle_ecl + atan(sinalpha, cosalpha)
	return (theta_gal, phi_gal, pol_angle_gal)
end

"""
    rotate_coordinates_e2g!(θ, φ, ψ)

Rotate arrays of pointings and orientation angles in place from Ecliptic to
Galactic coordinates.
"""
function rotate_coordinates_e2g!(θ, ϕ, ψ)
	@inbounds for i in eachindex(θ)
		θ[i], ϕ[i], ψ[i] = _rotate_coordinates_and_pol_e2g_for_one_sample(θ[i], ϕ[i], ψ[i])
	end
end
