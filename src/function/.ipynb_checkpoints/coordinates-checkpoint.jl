const EQUINOX = 2000

function thetaPhi2RaDec(θ, ϕ)
    dec =  π/2 - θ
    ra  = ϕ
    return (ra, dec)
end

function rotate_vec_ecl_to_gal_coordinates!(v::Vector{<:AbstractFloat},; equinox = EQUINOX)
    θ, ϕ     = vec2ang(v[1], v[2], v[3])
    ra, dec  = thetaPhi2RaDec(θ, ϕ)
    ecliptic = EclipticCoords{equinox}(ra, dec)
    gal      = convert(GalCoords, ecliptic)
    θ_gal, ϕ_gal = π/2 - gal.b, gal.l
    v_gal    = ang2vec(θ_gal, ϕ_gal)
    v       .= v_gal
    return v
end

x_axis_vec_gal = rotate_vec_ecl_to_gal_coordinates!([1.,0.,0.], equinox = EQUINOX) 
y_axis_vec_gal = rotate_vec_ecl_to_gal_coordinates!([0.,1.,0.], equinox = EQUINOX) 
z_axis_vec_gal = rotate_vec_ecl_to_gal_coordinates!([0.,0.,1.], equinox = EQUINOX)

const ECL_TO_GAL_ROT_MATRIX = SMatrix{3,3,Float64}([x_axis_vec_gal y_axis_vec_gal z_axis_vec_gal])
const NORTH_POLE_VEC = SVector{3,Float64}(z_axis_vec_gal)

function _ang2galvec_one_sample(theta, phi)
    """Transform a direction (theta, phi) in Ecliptic coordinates to
    a unit vector in Galactic coordinates.

    Parameters
    ----------
    theta : float, scalar
      The angle θ (colatitude) in Ecliptic coordinates
    phi : float, scalar
      The angle φ (longitude) in Ecliptic coordinates

    Returns
    -------
    vx, vy, vz : float
      A tuple of three floats

    See Also
    --------
    https://github.com/healpy/healpy/blob/main/healpy/rotator.py#L657
    """
    rotmatr = ECL_TO_GAL_ROT_MATRIX
    st = sin(theta)
    v =  SVector{3,Float64}([st * cos(phi); st * sin(phi); cos(theta)])
    return rotmatr * v
end
#=
function _vec2ang_for_one_sample(vx, vy, vz)
    """Transform a vector to angle given by theta,phi.

    Parameters
    ----------
    vx : float, scalar
      The x component of the vector (scalar)
    vy : float, scalar
      The y component of the vector (scalar))
    vz : float, scalar
      The z component of the vector (scalar)

    Returns
    -------
    theta, phi : float
      A tuple containing the value of the colatitude and of the longitude

    See Also
    --------
    https://github.com/healpy/healpy/blob/main/healpy/rotator.py#L610
    """

    return (atan(sqrt(vx^2 + vy^2), vz), atan(vy, vx))
end
=#
function _rotate_coordinates_and_pol_e2g_for_one_sample(theta_ecl, phi_ecl, pol_angle_ecl)
    """
    Rotate the angles theta,phi and psi from ecliptic to galactic coordinates

    Parameters
    ----------
    theta_ecl : scalar
      latitude (in radians) in the ecliptic coordinates
    phi_ecl   : scalar
      longitude (in radians) in the ecliptic coordinates
      
    pol_angle_ecl : scalar
      polarization angle (in radians) in ecliptic coordinates

    Returns
    -------
    theta_gal : scalar
      latitude (in radians) in the galactic coordinates
    phi_gal   : scalar
      longitude (in radians) in the galactic coordinates
      
    pol_angle_gal : scalar
      polarization angle (in radians) in galactic coordinates
    """

    vec = _ang2galvec_one_sample(theta_ecl, phi_ecl)
    theta_gal, phi_gal = vec2ang(vec[1], vec[2], vec[3])

    sinalpha = NORTH_POLE_VEC[1] * vec[2] - NORTH_POLE_VEC[2] * vec[1]
    cosalpha = NORTH_POLE_VEC[3] - vec[3] * (NORTH_POLE_VEC ⋅ vec)
    #* (
    #    NORTH_POLE_VEC[1] * vec[1]
    #    + NORTH_POLE_VEC[2] * vec[2]
    #    + NORTH_POLE_VEC[3] * vec[3]
    #)
    pol_angle_gal = pol_angle_ecl + atan(sinalpha, cosalpha)

    return (theta_gal, phi_gal, pol_angle_gal)
end

function rotate_coordinates_e2g!(θ, ϕ, ψ)
    @inbounds for i in eachindex(θ)
        θ[i], ϕ[i], ψ[i] = _rotate_coordinates_and_pol_e2g_for_one_sample(θ[i], ϕ[i], ψ[i])
    end
end
    