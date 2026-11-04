# Detector quaternions at the samples of an observation.
#
# Exact mode (delta_time_s === nothing), for every sample time t:
#     q(t) = G · ((A(t) · S(t)) · b) · d(t)
# Interpolated mode (delta_time_s = Δ), as litebird_sim:
#     B_j = G · ((A(t_j) · S(t_j)) · b)  on anchors t_j = t_sim + j Δ
#     q(t) = slerp(B · d)(t)               if d is constant
#     q(t) = slerp(B)(t) · slerp(d)(t)     otherwise
# with A the attitude, S the spin perturbation, b = bore2spin, d the detector
# quaternion and G the Ecliptic-to-output rotation.

struct _PointingPlan{O <: Observation}
	obs::O
	coord::CoordinateSystem
	g::Quaternion{Float64}
	b::Quaternion{Float64}
	detector_idx::Vector{Int}
	# interpolated mode only
	anchor_t0::Float64
	anchor_rate::Float64
	anchors::Vector{Vector{Quaternion{Float64}}}  # per selected detector (B·d, or B shared)
end

@inline _attitude_quat(a::SpinningScanningStrategy, t) = spin_to_ecliptic(a, t)
@inline _attitude_quat(a::RotQuaternion, t) = slerp(a, t)

@inline _perturbed(A::Quaternion, ::Nothing, _) = A
@inline _perturbed(A::Quaternion, s::RotQuaternion, t) = A * slerp(s, t)

@inline function _boresight_quat(obs::Observation, coord::CoordinateSystem, g, b, t)
	bore = _perturbed(_attitude_quat(obs.attitude, t), obs.spin_perturbation, t) * b
	return coord == Ecliptic ? bore : g * bore
end

function _detector_indices(obs::Observation, detector_idx)
	idx = detector_idx === :all ? collect(1:n_detectors(obs)) :
		  detector_idx isa Integer ? [Int(detector_idx)] : collect(Int, detector_idx)
	for i in idx
		1 <= i <= n_detectors(obs) || throw(ArgumentError("detector index $i is out of range 1:$(n_detectors(obs))"))
	end
	return idx
end

function _check_time_coverage(obs::Observation, t_first, t_last)
	if obs.attitude isa RotQuaternion
		check_covers(obs.attitude, t_first, t_last, "the attitude")
	end
	if obs.spin_perturbation !== nothing
		check_covers(obs.spin_perturbation, t_first, t_last, "the spin perturbation")
	end
	return nothing
end

function _pointing_plan(obs::Observation, coord::CoordinateSystem, detector_idx)
	idx = _detector_indices(obs, detector_idx)
	g = coordinate_quat(coord)
	b = bore2spin_quat(obs.instrument)
	t_first, t_last = sample_time(obs, 1), sample_time(obs, obs.n_samples)
	for i in idx
		check_covers(obs.detectors[i].quat, t_first, t_last, "the quaternion of detector \"$(obs.detectors[i].name)\"")
	end
	if obs.delta_time_s === nothing
		_check_time_coverage(obs, t_first, t_last)
		return _PointingPlan(obs, coord, g, b, idx, 0.0, 0.0, Vector{Quaternion{Float64}}[])
	end

	Δ = obs.delta_time_s
	tsim = obs.start_time_simulation_s
	j_lo = floor(Int, (t_first - tsim) / Δ)
	j_hi = max(ceil(Int, (t_last - tsim) / Δ), j_lo + 1)
	anchor_times = [tsim + j * Δ for j in j_lo:j_hi]
	_check_time_coverage(obs, anchor_times[1], anchor_times[end])
	bore = Vector{Quaternion{Float64}}(undef, length(anchor_times))
	@threads for j in eachindex(anchor_times)
		bore[j] = _boresight_quat(obs, coord, g, b, anchor_times[j])
	end
	anchors = map(idx) do i
		d = obs.detectors[i].quat
		is_constant(d) ? [normalize(q * d.quats[1]) for q in bore] : bore
	end
	return _PointingPlan(obs, coord, g, b, idx, anchor_times[1], 1 / Δ, anchors)
end

@inline function _interp_anchor(anchors, t0, rate, t)
	fi = (t - t0) * rate
	i = clamp(floor(Int, fi), 0, length(anchors) - 2)
	return slerp(anchors[i+1], anchors[i+2], fi - i)
end

"Quaternion of the `n`-th selected detector at sample `k`."
@inline function _sample_quat(p::_PointingPlan, k::Integer, n::Integer)
	obs = p.obs
	t = sample_time(obs, k)
	d = obs.detectors[p.detector_idx[n]].quat
	if obs.delta_time_s === nothing
		bore = _boresight_quat(obs, p.coord, p.g, p.b, t)
		return is_constant(d) ? bore * d.quats[1] : bore * slerp(d, t)
	end
	q = _interp_anchor(p.anchors[n], p.anchor_t0, p.anchor_rate, t)
	return is_constant(d) ? q : q * slerp(d, t)
end

@inline _hwp_angle_at(::Observation{<:Any, Nothing}, _) = 0.0
# Subtract the start times first: the absolute sample time loses ~1e-10 s of
# precision at t ~ 1e6 s, which ω t would amplify.
@inline _hwp_angle_at(obs::Observation{<:Any, IdealHWP}, k) =
	hwp_angle(obs.hwp, (obs.start_time_s - obs.start_time_simulation_s) + (k - 1) * (1 / obs.sampling_rate_hz))

_hwp_angles(::Observation{<:Any, Nothing}) = nothing
_hwp_angles(obs::Observation{<:Any, IdealHWP}) = [_hwp_angle_at(obs, k) for k in 1:obs.n_samples]

"""
    get_pointings(obs; detector_idx = :all, output_coordinate_system = Ecliptic) -> (pointings, hwp_angle)

Pointings of the detectors of `obs`, as `Observation.get_pointings` in
litebird_sim.

`pointings[:, k, i]` holds `(θ, φ, ψ)` (rad) of the `i`-th selected detector at
the `k`-th sample: colatitude, longitude in [0, 2π) and orientation angle in
(−π, π] (see [`compute_pointing_and_orientation`](@ref)). With
`detector_idx = :all` or a vector of (1-based) indices the array has size
`(3, n_samples, n_detectors)`; with an integer it has size `(3, n_samples)`.
The memory layout is the same as litebird_sim's `(n_detectors, n_samples, 3)`.

`hwp_angle` is the HWP angle of every sample, or `nothing` without an HWP.
"""
function get_pointings(obs::Observation; detector_idx = :all, output_coordinate_system::CoordinateSystem = Ecliptic)
	plan = _pointing_plan(obs, output_coordinate_system, detector_idx)
	ndet = length(plan.detector_idx)
	pt = Array{Float64}(undef, 3, obs.n_samples, ndet)
	@threads for k in 1:obs.n_samples
		@inbounds for n in 1:ndet
			θ, φ, ψ = compute_pointing_and_orientation(_sample_quat(plan, k, n))
			pt[1, k, n] = θ
			pt[2, k, n] = φ
			pt[3, k, n] = ψ
		end
	end
	out = detector_idx isa Integer ? pt[:, :, 1] : pt
	return out, _hwp_angles(obs)
end

"""
    PointingPixels

Result of [`get_pointing_pixels`](@ref): `pixels` (1-based HEALPix indices) and
`psi` (orientation angles, rad), both `(n_samples, n_detectors)`, the HWP angle
of every sample (`nothing` without an HWP), and the `nside`, `nest` ordering,
`coordinate_system` and `detectors` names they refer to.
"""
struct PointingPixels
	pixels::Matrix{Int}
	psi::Matrix{Float64}
	hwp_angle::Union{Nothing, Vector{Float64}}
	nside::Int
	nest::Bool
	coordinate_system::CoordinateSystem
	detectors::Vector{String}
end

@inline _ang2pix(resol::Resolution, nest::Bool, θ, φ) = nest ? ang2pixNest(resol, θ, φ) : ang2pixRing(resol, θ, φ)

"""
    get_pointing_pixels(obs; nside, nest = false, output_coordinate_system = Galactic,
                        detector_idx = :all) -> PointingPixels

HEALPix pixel indices (1-based, as in Healpix.jl; subtract 1 for healpy),
orientation angles and HWP angles of every sample, ready for map-making. The
pixels are computed sample by sample, without storing `(θ, φ)`.
"""
function get_pointing_pixels(obs::Observation; nside::Integer, nest::Bool = false,
	output_coordinate_system::CoordinateSystem = Galactic, detector_idx = :all)
	plan = _pointing_plan(obs, output_coordinate_system, detector_idx)
	resol = Resolution(nside)
	ndet = length(plan.detector_idx)
	pixels = Matrix{Int}(undef, obs.n_samples, ndet)
	psi = Matrix{Float64}(undef, obs.n_samples, ndet)
	@threads for k in 1:obs.n_samples
		@inbounds for n in 1:ndet
			θ, φ, ψ = compute_pointing_and_orientation(_sample_quat(plan, k, n))
			pixels[k, n] = _ang2pix(resol, nest, θ, φ)
			psi[k, n] = ψ
		end
	end
	names = [obs.detectors[i].name for i in plan.detector_idx]
	return PointingPixels(pixels, psi, _hwp_angles(obs), nside, nest, output_coordinate_system, names)
end

const _CHUNK_SAMPLES = 1 << 18

"""
    _foreach_chunk(f, plan, resol, nest)

Call `f(first_sample, pixels, psi)` for consecutive chunks of samples, where
`pixels` and `psi` are `(chunk_length, n_selected_detectors)` buffers.
"""
function _foreach_chunk(f::F, plan::_PointingPlan, resol::Resolution, nest::Bool) where {F}
	obs = plan.obs
	ndet = length(plan.detector_idx)
	len = min(_CHUNK_SAMPLES, obs.n_samples)
	pixels = Matrix{Int}(undef, len, ndet)
	psi = Matrix{Float64}(undef, len, ndet)
	for first in 1:len:obs.n_samples
		m = min(len, obs.n_samples - first + 1)
		@threads for j in 1:m
			@inbounds for n in 1:ndet
				θ, φ, ψ = compute_pointing_and_orientation(_sample_quat(plan, first + j - 1, n))
				pixels[j, n] = _ang2pix(resol, nest, θ, φ)
				psi[j, n] = ψ
			end
		end
		f(first, view(pixels, 1:m, :), view(psi, 1:m, :))
	end
	return nothing
end

"""
    get_hitmap(obs; nside, output_coordinate_system = Galactic, nest = false,
               detector_idx = :all) -> Vector{Int}

Number of samples falling in every HEALPix pixel, summed over the selected
detectors (and over the observations if `obs` is a vector), as `get_hitmap` in
litebird_sim.
"""
function get_hitmap(obs::Observation; nside::Integer, output_coordinate_system::CoordinateSystem = Galactic,
	nest::Bool = false, detector_idx = :all)
	return get_hitmap([obs]; nside, output_coordinate_system, nest, detector_idx)
end

function get_hitmap(observations::AbstractVector{<:Observation}; nside::Integer,
	output_coordinate_system::CoordinateSystem = Galactic, nest::Bool = false, detector_idx = :all)
	resol = Resolution(nside)
	hits = zeros(Int, resol.numOfPixels)
	for obs in observations
		plan = _pointing_plan(obs, output_coordinate_system, detector_idx)
		_foreach_chunk(plan, resol, nest) do _, pixels, _
			@inbounds for p in pixels
				hits[p] += 1
			end
		end
	end
	return hits
end
