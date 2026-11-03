const EX = SVector(1.0, 0.0, 0.0)
const EY = SVector(0.0, 1.0, 0.0)
const EZ = SVector(0.0, 0.0, 1.0)

"Angular velocity of the anti-Sun direction (rad/s), one revolution per 365.25 days."
const OMEGA_REVOL = 2π / (60.0 * 60.0 * 24.0 * 365.25)

"""
    rotate_quat(omega, t, axis) -> Quaternion

Quaternion of the rotation by `mod2pi(omega * t)` around `axis`.
"""
@inline rotate_quat(omega, t, rotate_axis) = rotation_quat(rotate_axis, mod2pi(omega * t))

"""
    rotate_quat(ϕ, axis) -> Quaternion

Quaternion of the rotation by `ϕ` (rad) around `axis`.
"""
@inline rotate_quat(ϕ, rotate_axis) = rotation_quat(rotate_axis, ϕ)

mutable struct ScanningStrategy{T <: AbstractFloat, I <: Int, AS <: AbstractString}
	nside::I
	duration::I
	sampling_rate::T
	alpha::T
	beta::T
	gamma::T
	prec_rpm::T
	spin_rpm::T
	hwp_rpm::T
	start_point::AS
	start_angle::T
	coord::AS
	quat::Vector{Vector{Float64}}
	name::Vector{String}
	info::DataFrame
end

function show_ss(ss::ScanningStrategy)
	@printf("%-24s : %i \n", "nside", ss.nside)
	@printf("%-24s : %0.1f \n", "duration [sec]", ss.duration)
	@printf("%-24s : %0.1f \n", "sampling rate [Hz]", ss.sampling_rate)
	@printf("%-24s : %0.1f \n", "alpha [deg]", ss.alpha)
	@printf("%-24s : %0.1f \n", "beta [deg]", ss.beta)
	@printf("%-24s : %0.1f \n", "gamma [deg]", ss.gamma)

	@printf("%-24s : %0.3f\n", "prec. period [min]", period2rpm(ss.prec_rpm))
	@printf("%1s %-22s : %f\n", '↳', "prec. rate [rpm]", ss.prec_rpm)

	@printf("%-24s : %0.3f\n", "spin period [min]", period2rpm(ss.spin_rpm))
	@printf("%1s %-22s : %f\n", '↳', "spin rate [rpm]", ss.spin_rpm)
	@printf("%-24s : %f \n", "HWP rot. rate[rpm]", ss.hwp_rpm)
	@printf("%-24s : %f \n", "start angle", ss.start_angle)
	@printf("%-24s : %s \n", "coordinate system", ss.coord)
	@printf("%-24s\n", "FPU")
	for i in eachindex(ss.quat)
		@printf("↳ Det. %i  name | %s %-12s : (x,y,z,w) = [%0.3f, %0.3f, %0.3f, %0.3f] \n", i, ss.name[i], "",
			ss.quat[i][1],
			ss.quat[i][2],
			ss.quat[i][3],
			ss.quat[i][4]
		)
	end
end

"""
    Pointings{T}

Pointing information of one sample: unit vector `(x, y, z)`, sky position
`(θ, φ)`, RING pixel index `Ω`, orientation angle `ψ`, HWP angle `ϕ` and
`ξ = mod2pi(2ϕ) + ψ`.
"""
struct Pointings{T <: AbstractFloat}
	x::T
	y::T
	z::T
	θ::T
	φ::T
	Ω::Int
	ψ::T
	ϕ::T
	ξ::T
end

function pointings(resol::Resolution, θ, φ, ψ, ϕ)
	x, y, z = ang2vec(θ, φ)
	Ω = ang2pixRing(resol, θ, φ)
	ξ = mod2pi(2ϕ) + ψ
	T = float(promote_type(typeof(θ), typeof(φ), typeof(ψ), typeof(ϕ)))
	return Pointings{T}(x, y, z, θ, φ, Ω, ψ, ϕ, ξ)
end

function gen_ScanningStrategy(;
	nside = 128,
	duration = 60 * 60 * 24 * 365,
	sampling_rate = 1.0,
	alpha = 45,
	beta = 50,
	gamma = 0,
	prec_rpm = period2rpm(192.348),
	spin_rpm = 0.05,
	hwp_rpm = 0,
	start_point = "equator",
	start_angle = 0.0,
	coord = "E",
	quat = [[0.0, 0.0, 1.0, 0.0]],
	name = ["boresight"],
	info = DataFrame(),
)
	ScanningStrategy(
		nside,
		duration,
		sampling_rate,
		Float64(alpha),
		Float64(beta),
		Float64(gamma),
		Float64(prec_rpm),
		Float64(spin_rpm),
		Float64(hwp_rpm),
		start_point,
		start_angle,
		coord,
		quat,
		name,
		info,
	)
end

function get_satellite(satellite)
	six_month_min = 24 * 60 * 60 * (365.25 / 2) / 60
	four_day_min  = 24 * 60 * 60 * 4 / 60
	satellites    = ["WMAP", "Planck", "PICO", "CORE", "EPIC", "LiteBIRD"]
	Alpha         = [70.0, 7.5, 26, 30, 45, 45]
	Beta          = [22.5, 85, 69, 65, 55, 50]
	T_alpha       = [60.0, six_month_min, 10 * 60, four_day_min, 3.2 * 60, 192.348] # min
	T_beta        = [129 / 60, 1, 1, 2, 1, 20.0] # min

	if satellite in satellites
		idx = findfirst(x -> x == satellite, satellites)
	else
		return "Satellite not found in the list."
	end
	ss = gen_ScanningStrategy(;
		nside = 128,
		duration = 60 * 60 * 24 * 365,
		sampling_rate = 1.0,
		alpha = Alpha[idx],
		beta = Beta[idx],
		gamma = 0,
		prec_rpm = period2rpm(T_alpha[idx]),
		spin_rpm = period2rpm(T_beta[idx]),
		hwp_rpm = 0,
		start_point = "equator",
		start_angle = 0.0,
		coord = "E",
		quat = [[0.0, 0.0, 1.0, 0.0]],
		name = ["boresight"],
		info = DataFrame(),
	)
end

_flip_angle(ss::ScanningStrategy) = ss.start_point == "pole" ? π : 0.0

"""
    imo2ecl_coordinates(ss) -> (pointing_dirs, pol_dirs)

Pointing and polarization directions of every detector in the spin-axis frame
at `t = 0`, as vectors of `SVector{3}`. Detectors that match none of the
supported naming schemes (boresight, IMO entries with a `pixtype` column, or
"999"-prefixed names) are left as zero vectors.
"""
function imo2ecl_coordinates(ss::ScanningStrategy)
	q_dets     = [zero(SVector{3, Float64}) for i in eachindex(ss.quat)]
	q_pol_dets = [zero(SVector{3, Float64}) for i in eachindex(ss.quat)]
	for i in eachindex(ss.quat)
		if ss.name[i] == "boresight"
			Q = rotate_quat(deg2rad(90.0 - ss.alpha), EY) * rotate_quat(_flip_angle(ss), EZ) * rotate_quat(deg2rad(ss.beta), EY)
		elseif hasproperty(ss.info, :pixtype) # IMO before the reformation
			telescope = split(ss.name[i], "_")[1]
			q_imo     = Quaternion(ss.quat[i])
			q_gamma   = rotate_quat(deg2rad(ss.gamma), EZ)
			Q = rotate_quat(deg2rad(ss.beta), EY) * q_gamma * q_imo
			if telescope == "001" || telescope == "002" # MFT, HFT
				Q = rotate_quat(π, EZ) * Q
			end
			Q = rotate_quat(deg2rad(90.0 - ss.alpha), EY) * Q
		elseif split(ss.name[i], "_")[1] == "999" # IMO after the reformation
			q_imo   = Quaternion(ss.quat[i])
			q_gamma = rotate_quat(deg2rad(90 + ss.gamma), EZ)
			Q = rotate_quat(deg2rad(90.0 - ss.alpha), EY) * (rotate_quat(deg2rad(ss.beta), EY) * q_gamma * q_imo)
		else
			continue
		end
		q_dets[i]     = rotate(Q, EZ)
		q_pol_dets[i] = rotate(Q, EX)
	end
	return (q_dets, q_pol_dets)
end

function _clip_sincos(x)
	# np.clip would be perfect here, but it is slightly inefficient to use here
	# because `x` might either be a float64 or a float32, and we would need to
	# define the −1 and 1 constant accordingly. Better to resort to a combination
	# of calls to min/max, as they are efficiently optimized by LLVM.
	return min(max(x, -1), 1)
end

function polarization_angle(θ, ϕ, poldir)
	cos_psi = _clip_sincos(-sin(ϕ) * poldir[1] + cos(ϕ) * poldir[2])
	sin_psi = _clip_sincos(
		(-cos(θ) * cos(ϕ) * poldir[1])
		+ (-cos(θ) * sin(ϕ) * poldir[2])
		+ (sin(θ) * poldir[3]),
	)
	return atan(sin_psi, cos_psi)
end

function _time_array(ss::ScanningStrategy, start, stop)
	if start > stop - 1 / ss.sampling_rate
		error("ERROR: \n The `start` time of the calculation is greater than or equal to the `stop` time.")
	end
	return start:1/ss.sampling_rate:stop-1/ss.sampling_rate |> LinRange
end

"""
    _pointings(ss, time_array, det_dirs, pol_dirs, revol_angle)

Compute (θ, φ, ψ) for every sample and detector. `revol_angle(t)` gives the
rotation of the anti-Sun direction at time `t` on top of `ss.start_angle`.
"""
function _pointings(ss::ScanningStrategy, time_array, det_dirs, pol_dirs, revol_angle::F) where {F}
	omega_spin = rpm2angfreq(ss.spin_rpm)
	omega_prec = rpm2angfreq(ss.prec_rpm)
	spin_axis  = SVector(cosd(ss.alpha), 0.0, sind(ss.alpha))
	q_sun      = rotate_quat(ss.start_angle, EZ)

	nsamples  = length(time_array)
	theta_tod = zeros(nsamples, length(det_dirs))
	phi_tod   = zeros(nsamples, length(det_dirs))
	psi_tod   = zeros(nsamples, length(det_dirs))
	@inbounds for i in eachindex(det_dirs)
		d, p = det_dirs[i], pol_dirs[i]
		@threads for j in 1:nsamples
			t = time_array[j]
			Q = q_sun * rotate_quat(revol_angle(t), EZ) * rotate_quat(omega_prec, t, EX) * rotate_quat(omega_spin, t, spin_axis)
			vec_point = rotate(Q, d)
			θ, ϕ = vec2ang(vec_point[1], vec_point[2], vec_point[3])
			theta_tod[j, i] = θ
			phi_tod[j, i] = ϕ
			psi_tod[j, i] = polarization_angle(θ, ϕ, rotate(Q, p))
		end
	end
	if ss.coord == "G"
		rotate_coordinates_e2g!(theta_tod, phi_tod, psi_tod)
	end
	return (theta_tod, phi_tod, psi_tod, time_array)
end

"""
    get_pointings(ss, start, stop) -> (θ, φ, ψ, time)

Pointing angles of every detector for `start ≤ t < stop` (s), sampled at
`ss.sampling_rate`. `θ`, `φ` and `ψ` have shape `(nsamples, ndetectors)`.
"""
function get_pointings(ss::ScanningStrategy, start, stop)
	time_array = _time_array(ss, start, stop)
	det_dirs, pol_dirs = imo2ecl_coordinates(ss)
	return _pointings(ss, time_array, det_dirs, pol_dirs, t -> mod2pi(OMEGA_REVOL * t))
end

"""
    get_pointings(ss, start, stop, step) -> (θ, φ, ψ, time)

Same as `get_pointings(ss, start, stop)`, but the anti-Sun direction advances
by 360/365.25 degrees every `step` seconds instead of continuously.
"""
function get_pointings(ss::ScanningStrategy, start, stop, step)
	time_array = _time_array(ss, start, stop)
	det_dirs, pol_dirs = imo2ecl_coordinates(ss)
	deg_per_day = 360 / 365.25
	return _pointings(ss, time_array, det_dirs, pol_dirs, t -> deg2rad(deg_per_day * floor(t / step)))
end

"""
	get_pointing_pixels(ss::ScanningStrategy, start, stop)

This function will return pointing pixel tod as tuple.
	# Arguments
	...
	- `ss::ScanningStrategy`: the ScanningStrategy struct.
	- `start::Int`: the initial time for pointing calculation.
	- `stop::Int`: the finish time for pointing calculation.
	...

	# Examples
	```jldoctest
	julia> ss = gen_ScanningStrategy()
	julia> pointings = get_pointing_pixels(s, 0, 100)

	# Returns
	...
	- `pointings[1]`: TOD of pixel indicies. shape:(duration*sampling_rate, numOfdet)
	- `pointings[2]`: TOD of psi. shape:(duration*sampling_rate, numOfdet)
	- `pointings[3]`: Array of time. shape:(duration*sampling_rate)
	...
"""
@inline function get_pointing_pixels(ss, start, stop)
	theta_tod, phi_tod, psi_tod, time_array = get_pointings(ss, start, stop)
	resol = Resolution(ss.nside)
	pix_tod = zeros(Int64, size(theta_tod))
	for i in eachindex(theta_tod)
		pix_tod[i] = ang2pixRing(resol, theta_tod[i], phi_tod[i])
	end
	return (pix_tod, psi_tod, time_array)
end

function angtod2hitmap(nside::Int, theta_tod, phi_tod)
	resol = Resolution(nside)
	hit_map = zeros(resol.numOfPixels)
	@inbounds for k in eachindex(theta_tod)
		hit_map[ang2pixRing(resol, theta_tod[k], phi_tod[k])] += 1
	end
	return hit_map
end

function pixtod2hitmap(nside::Int, pixtod)
	resol = Resolution(nside)
	hit_map = zeros(resol.numOfPixels)
	for i in eachindex(pixtod)
		hit_map[pixtod[i]] += 1
	end
	return hit_map
end

function normarize!(resol::Resolution, maps::Array, hitmap::Array)
	if size(maps) == (3, 1, resol.numOfPixels)
		for i in 1:3
			maps[i, 1, :] ./= hitmap
		end
	end
	if size(maps) == (3, 3, resol.numOfPixels)
		for i in 1:3
			for j in 1:3
				maps[i, j, :] ./= hitmap
			end
		end
	end
	return maps
end


# Pointings of detectors placed at arbitrary positions on the focal plane.

"""
    ang_to_quat(theta, phi) -> Vector{Float64}

Quaternion `(x, y, z, w)` that moves the boresight to the focal-plane position
`(theta, phi)` (rad).
"""
function ang_to_quat(theta, phi)
	q = rotate_quat(phi, EZ) * rotate_quat(theta, EY)
	return [q.x, q.y, q.z, q.w]
end

function get_q_dets_q_pol_dets(ss::ScanningStrategy)
	q_dets     = [zero(SVector{3, Float64}) for i in eachindex(ss.quat)]
	q_pol_dets = [zero(SVector{3, Float64}) for i in eachindex(ss.quat)]
	for i in eachindex(ss.quat)
		q_initial = Quaternion(ss.quat[i])
		Q = rotate_quat(deg2rad(90.0 - ss.alpha), EY) * rotate_quat(_flip_angle(ss), EZ) * rotate_quat(deg2rad(ss.beta), EY) * q_initial
		q_dets[i]     = rotate(Q, EZ)
		q_pol_dets[i] = rotate(Q, EX)
	end
	return (q_dets, q_pol_dets)
end

"""
    get_pointings_arbitrary_dets(ss, start, stop) -> (θ, φ, ψ, time)

Like `get_pointings`, but every entry of `ss.quat` is used directly as the
detector's focal-plane rotation (see [`ang_to_quat`](@ref)), whatever its name.
"""
function get_pointings_arbitrary_dets(ss::ScanningStrategy, start, stop)
	time_array = _time_array(ss, start, stop)
	det_dirs, pol_dirs = get_q_dets_q_pol_dets(ss)
	return _pointings(ss, time_array, det_dirs, pol_dirs, t -> mod2pi(OMEGA_REVOL * t))
end
