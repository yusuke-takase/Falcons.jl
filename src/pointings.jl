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
#
# Everything a sample needs is packed into a `_PlanData`, whose arrays are
# `Vector`s on the CPU and device arrays in the CUDA extension, so that both
# backends evaluate the same `_sample_quat`. Time coverage is checked on the
# host before a plan is built, so the evaluation itself never throws. In
# interpolated mode the anchors are computed on the host for one block of
# samples at a time (`_anchor_block`), as litebird_sim slices its quaternions
# per block, so their memory does not grow with the observation.

"A time series of quaternions sampled at `t0 + (k - 1) / rate`."
struct _Series{V}
	quats::V
	t0::Float64
	rate::Float64
end

_Series(r::RotQuaternion) = _Series(r.quats, r.start_time_s, r.sampling_rate_hz)

@inline function _interp(quats, len, offset, t0, rate, t)
	fi = (t - t0) * rate
	i = clamp(_ifloor(fi), 0, len - 2)
	return slerp(quats[offset+i+1], quats[offset+i+2], fi - i)
end

@inline _series_at(s::_Series, t) = _interp(s.quats, length(s.quats), 0, s.t0, s.rate, t)

@inline _attitude_at(a::SpinningScanningStrategy, t) = spin_to_ecliptic(a, t)
@inline _attitude_at(a::_Series, t) = _series_at(a, t)

@inline _perturbed(A::Quaternion, ::Nothing, _) = A
@inline _perturbed(A::Quaternion, s::_Series, t) = A * _series_at(s, t)

"Detector quaternions: constant ones, and time series concatenated in `series`."
struct _DetectorQuats{V, I, F}
	constant::V
	series::V
	offset::I
	len::I    # 0 for a constant quaternion
	t0::F
	rate::F
end

@inline _is_constant_det(d::_DetectorQuats, n) = d.len[n] == 0

@inline function _det_quat(d::_DetectorQuats, n, t)
	len = d.len[n]
	len == 0 && return d.constant[n]
	return _interp(d.series, len, d.offset[n], d.t0[n], d.rate[n], t)
end

function _DetectorQuats(rots::Vector{<:RotQuaternion})
	constant = Quaternion{Float64}[]
	series = Quaternion{Float64}[]
	offset, len, t0, rate = Int[], Int[], Float64[], Float64[]
	for r in rots
		push!(offset, length(series))
		if is_constant(r)
			push!(constant, r.quats[1])
			push!(len, 0); push!(t0, 0.0); push!(rate, 0.0)
		else
			push!(constant, one(Quaternion{Float64}))
			append!(series, r.quats)
			push!(len, length(r.quats)); push!(t0, r.start_time_s); push!(rate, r.sampling_rate_hz)
		end
	end
	isempty(series) && push!(series, one(Quaternion{Float64}))  # keep device arrays non-empty
	return _DetectorQuats(constant, series, offset, len, t0, rate)
end

struct _PlanData{A, P, D, M}
	attitude::A          # SpinningScanningStrategy or _Series
	perturbation::P      # Nothing or _Series
	b::Quaternion{Float64}
	g::Quaternion{Float64}
	galactic::Bool
	dets::D
	interp::Bool
	# Interpolated mode: anchor i = 0, …, n_anchors − 1 of the observation is at
	# anchor_tsim + (anchor_j0 + i) anchor_dt = anchor_t0 + i / anchor_rate, and
	# anchors[:, n] holds anchors anchor_first, anchor_first + 1, … of detector n.
	anchors::M
	anchor_first::Int
	n_anchors::Int
	anchor_j0::Int
	anchor_tsim::Float64
	anchor_dt::Float64
	anchor_t0::Float64
	anchor_rate::Float64
	start_time_s::Float64
	dt::Float64          # 1 / sampling_rate_hz
	has_hwp::Bool
	hwp_omega::Float64
	hwp_start::Float64
	hwp_dt0::Float64     # start_time_s − start_time_simulation_s
end

# Products followed by sums are written as explicit fma: the GPU compiler fuses
# them anyway, so this keeps the CPU and GPU results bit-identical.
@inline _plan_time(p::_PlanData, k) = fma(Float64(k - 1), p.dt, p.start_time_s)

@inline function _boresight_at(p::_PlanData, t)
	bore = _perturbed(_attitude_at(p.attitude, t), p.perturbation, t) * p.b
	return p.galactic ? p.g * bore : bore
end

"Index (0-based) of the anchor at or before time `t`, and the position of `t` from anchor 0."
@inline function _anchor_index(p::_PlanData, t)
	fi = (t - p.anchor_t0) * p.anchor_rate
	return clamp(_ifloor(fi), 0, p.n_anchors - 2), fi
end

"Quaternion of the `n`-th selected detector at sample `k` (in interpolated mode, `p` must hold its anchors)."
@inline function _sample_quat(p::_PlanData, k, n)
	t = _plan_time(p, k)
	if !p.interp
		return _boresight_at(p, t) * _det_quat(p.dets, n, t)
	end
	i, fi = _anchor_index(p, t)
	j = i - p.anchor_first
	q = slerp(p.anchors[j+1, n], p.anchors[j+2, n], fi - i)
	return _is_constant_det(p.dets, n) ? q : q * _det_quat(p.dets, n, t)
end

# Subtract the start times first: the absolute sample time loses ~1e-10 s of
# precision at t ~ 1e6 s, which ω t would amplify.
@inline _hwp_at(p::_PlanData, k) = p.has_hwp ? _rem2pi(fma(p.hwp_omega, fma(Float64(k - 1), p.dt, p.hwp_dt0), p.hwp_start)) : 0.0

function _detector_indices(obs::Observation, detector_idx)
	idx = detector_idx === :all ? collect(1:n_detectors(obs)) :
		  detector_idx isa Integer ? [Int(detector_idx)] : collect(Int, detector_idx)
	for i in idx
		1 <= i <= n_detectors(obs) || throw(ArgumentError("detector index $i is out of range 1:$(n_detectors(obs))"))
	end
	return idx
end

function _check_time_coverage(obs::Observation, t_first, t_last)
	obs.attitude isa RotQuaternion && check_covers(obs.attitude, t_first, t_last, "the attitude")
	obs.spin_perturbation !== nothing && check_covers(obs.spin_perturbation, t_first, t_last, "the spin perturbation")
	return nothing
end

_series_or_analytic(a::SpinningScanningStrategy) = a
_series_or_analytic(a::RotQuaternion) = is_constant(a) ? _Series(fill(a.quats[1], 2), 0.0, 0.0) : _Series(a)
_series_or_nothing(::Nothing) = nothing
_series_or_nothing(a::RotQuaternion) = _series_or_analytic(a)

"Throw an `ArgumentError` unless every quaternion of `quats` is a finite unit quaternion."
function _check_unit_quats(quats, what::AbstractString)
	for (k, q) in enumerate(quats)
		n2 = abs2(q)
		(isfinite(n2) && abs(n2 - 1) <= 1e-6) ||
			throw(ArgumentError("$what: quaternion $k is $q, not a finite unit quaternion (missing data?)"))
	end
	return nothing
end

function _check_finite(x::Real, what::AbstractString)
	isfinite(x) || throw(ArgumentError("$what is $x"))
	return nothing
end

"""
    _check_inputs(obs, idx)

Throw an `ArgumentError` for a non-finite time, rate or angle, or a zero,
non-finite or non-normalised quaternion (e.g. a missing sample of an
attitude): the GPU kernels cannot raise, and would otherwise accumulate
into undefined pixels.
"""
function _check_inputs(obs::Observation, idx)
	_check_finite(obs.start_time_s, "start_time_s")
	_check_finite(obs.start_time_simulation_s, "start_time_simulation_s")
	_check_finite(obs.sampling_rate_hz, "sampling_rate_hz")
	obs.delta_time_s === nothing || _check_finite(obs.delta_time_s, "delta_time_s")
	if obs.attitude isa SpinningScanningStrategy
		for f in fieldnames(SpinningScanningStrategy)
			_check_finite(getfield(obs.attitude, f), "the scanning strategy's $f")
		end
	else
		_check_unit_quats(obs.attitude.quats, "the attitude")
	end
	obs.spin_perturbation === nothing || _check_unit_quats(obs.spin_perturbation.quats, "the spin perturbation")
	_check_unit_quats((bore2spin_quat(obs.instrument),), "the boresight-to-spin rotation of the instrument")
	for i in idx
		_check_unit_quats(obs.detectors[i].quat.quats, "the quaternion of detector \"$(obs.detectors[i].name)\"")
	end
	if obs.hwp !== nothing
		_check_finite(obs.hwp.ang_speed_radpsec, "the HWP's ang_speed_radpsec")
		_check_finite(obs.hwp.start_angle_rad, "the HWP's start_angle_rad")
	end
	return nothing
end

const _NONFINITE_ANGLES = "some samples have non-finite pointing or HWP angles; check the quaternions, " *
						  "the scan parameters and the times"

"""
    _plan_data(obs, coord, detector_idx) -> (_PlanData, selected detector indices)

Check the inputs (see `_check_inputs`) and that every time series covers the
samples (or anchors) it is evaluated at, and precompute what the samples need.
"""
function _plan_data(obs::Observation, coord::CoordinateSystem, detector_idx)
	idx = _detector_indices(obs, detector_idx)
	_check_inputs(obs, idx)
	g = coordinate_quat(coord)
	b = bore2spin_quat(obs.instrument)
	t_first, t_last = sample_time(obs, 1), sample_time(obs, obs.n_samples)
	for i in idx
		check_covers(obs.detectors[i].quat, t_first, t_last, "the quaternion of detector \"$(obs.detectors[i].name)\"")
	end
	dets = _DetectorQuats([obs.detectors[i].quat for i in idx])
	attitude = _series_or_analytic(obs.attitude)
	perturbation = _series_or_nothing(obs.spin_perturbation)
	hwp = obs.hwp
	hwp_args = hwp === nothing ? (false, 0.0, 0.0, 0.0) :
			   (true, hwp.ang_speed_radpsec, hwp.start_angle_rad, obs.start_time_s - obs.start_time_simulation_s)
	common = (b, g, coord == Galactic, dets)
	no_anchors = Matrix{Quaternion{Float64}}(undef, 0, length(idx))
	times = (obs.start_time_s, 1 / obs.sampling_rate_hz)

	if obs.delta_time_s === nothing
		_check_time_coverage(obs, t_first, t_last)
		data = _PlanData(attitude, perturbation, common..., false, no_anchors, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, times...,
			hwp_args...)
		return data, idx
	end

	Δ = obs.delta_time_s
	tsim = obs.start_time_simulation_s
	j_lo = floor(Int, (t_first - tsim) / Δ)
	j_hi = max(ceil(Int, (t_last - tsim) / Δ), j_lo + 1)
	_check_time_coverage(obs, tsim + j_lo * Δ, tsim + j_hi * Δ)
	data = _PlanData(attitude, perturbation, common..., true, no_anchors, 0, j_hi - j_lo + 1, j_lo, tsim, Δ,
		tsim + j_lo * Δ, 1 / Δ, times..., hwp_args...)
	return data, idx
end

"A copy of `p` holding `anchors`, the first of which is anchor `first` (0-based)."
function _with_anchors(p::_PlanData, anchors, first::Int)
	return _PlanData(p.attitude, p.perturbation, p.b, p.g, p.galactic, p.dets, p.interp, anchors, first, p.n_anchors,
		p.anchor_j0, p.anchor_tsim, p.anchor_dt, p.anchor_t0, p.anchor_rate, p.start_time_s, p.dt, p.has_hwp,
		p.hwp_omega, p.hwp_start, p.hwp_dt0)
end

"""
    _anchor_block(p, k_first, k_last) -> _PlanData

A copy of the host plan `p` with only the anchors that the samples
`k_first:k_last` interpolate between (`B · d` for a constant detector
quaternion, `B` otherwise). The anchors are the same as if all of them had
been computed at once. `p` itself is returned in exact mode.
"""
function _anchor_block(p::_PlanData, k_first::Integer, k_last::Integer)
	p.interp || return p
	i_lo = first(_anchor_index(p, _plan_time(p, k_first)))
	i_hi = first(_anchor_index(p, _plan_time(p, k_last))) + 1
	bore = Vector{Quaternion{Float64}}(undef, i_hi - i_lo + 1)
	@threads for i in i_lo:i_hi
		bore[i-i_lo+1] = _boresight_at(p, p.anchor_tsim + (p.anchor_j0 + i) * p.anchor_dt)
	end
	ndet = length(p.dets.len)
	anchors = Matrix{Quaternion{Float64}}(undef, length(bore), ndet)
	for n in 1:ndet
		if _is_constant_det(p.dets, n)
			anchors[:, n] .= normalize.(bore .* Ref(p.dets.constant[n]))
		else
			anchors[:, n] .= bore
		end
	end
	return _with_anchors(p, anchors, i_lo)
end

# Samples per block on the CPU: bounds the anchors and buffers held at once.
const _CHUNK_SAMPLES = Ref(1 << 18)

"Samples per block: at most `len`, and about as many anchors."
function _block_length(p::_PlanData, nsamples::Integer, len::Integer = _CHUNK_SAMPLES[])
	if p.interp
		# Samples sparser than the anchors: bound the anchors of a block as well.
		len = min(len, max(1, floor(Int, len * p.anchor_dt / p.dt)))
	end
	return max(1, min(len, nsamples))
end

"""
    _foreach_block(f, p, nsamples)

Call `f(k_first, k_last, block_plan)` for consecutive blocks of samples, where
`block_plan` holds the anchors of the block (see `_anchor_block`).
"""
function _foreach_block(f::F, p::_PlanData, nsamples::Integer) where {F}
	len = _block_length(p, nsamples)
	for k_first in 1:len:nsamples
		k_last = min(k_first + len - 1, nsamples)
		f(k_first, k_last, _anchor_block(p, k_first, k_last))
	end
	return nothing
end

_hwp_angles(obs::Observation, p::_PlanData) = obs.hwp === nothing ? nothing : [_hwp_at(p, k) for k in 1:obs.n_samples]

# ----------------------------------------------------------------------------
# Backends

"""
    AbstractBackend

Where pointings are computed: [`CPUBackend`](@ref) or [`GPUBackend`](@ref).
"""
abstract type AbstractBackend end

"Compute on the CPU with `Threads.nthreads()` threads (the default)."
struct CPUBackend <: AbstractBackend end

"""
    GPUBackend(; num_gpus = 1, gpu_ids = nothing)

Compute on CUDA GPUs in Float64. Requires `using CUDA` (which loads the
`FalconsCUDAExt` extension). `gpu_ids` selects devices by their 0-based ids;
otherwise the first `num_gpus` devices are used, splitting the samples
between them. Start Julia with at least as many threads as GPUs so that the
devices run concurrently.
"""
struct GPUBackend <: AbstractBackend
	num_gpus::Int
	gpu_ids::Union{Nothing, Vector{Int}}
end

function GPUBackend(; num_gpus::Integer = 1, gpu_ids = nothing)
	num_gpus >= 1 || throw(ArgumentError("num_gpus must be ≥ 1, got $num_gpus"))
	ids = gpu_ids === nothing ? nothing : collect(Int, gpu_ids)
	ids === nothing || (!isempty(ids) && allunique(ids)) || throw(ArgumentError("gpu_ids must be non-empty and unique"))
	return GPUBackend(ids === nothing ? num_gpus : length(ids), ids)
end

# Implemented by the CUDA extension.
function _gpu_pointings end
function _gpu_pointing_pixels end
function _gpu_hits! end
function _gpu_h_maps! end

function _require_gpu()
	Base.get_extension(@__MODULE__, :FalconsCUDAExt) === nothing &&
		throw(ArgumentError("GPUBackend needs the CUDA extension: run `using CUDA` before calling Falcons with a GPU backend"))
	return nothing
end

# ----------------------------------------------------------------------------
# Public functions

"""
    get_pointings(obs; detector_idx = :all, output_coordinate_system = Ecliptic,
                  backend = CPUBackend()) -> (pointings, hwp_angle)

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
function get_pointings(obs::Observation; detector_idx = :all, output_coordinate_system::CoordinateSystem = Ecliptic,
	backend::AbstractBackend = CPUBackend())
	plan, idx = _plan_data(obs, output_coordinate_system, detector_idx)
	pt = if backend isa GPUBackend
		_require_gpu()
		_gpu_pointings(plan, obs.n_samples, length(idx), backend)
	else
		_cpu_pointings(plan, obs.n_samples, length(idx))
	end
	out = detector_idx isa Integer ? pt[:, :, 1] : pt
	return out, _hwp_angles(obs, plan)
end

function _cpu_pointings(plan::_PlanData, nsamples::Int, ndet::Int)
	pt = Array{Float64}(undef, 3, nsamples, ndet)
	_foreach_block(plan, nsamples) do k_first, k_last, bp
		@threads for k in k_first:k_last
			@inbounds for n in 1:ndet
				θ, φ, ψ = compute_pointing_and_orientation(_sample_quat(bp, k, n))
				pt[1, k, n] = θ
				pt[2, k, n] = φ
				pt[3, k, n] = ψ
			end
		end
	end
	return pt
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

@inline _healpix_ang2pix(resol::Resolution, nest::Bool, θ, φ) = nest ? ang2pixNest(resol, θ, φ) : ang2pixRing(resol, θ, φ)

"""
    get_pointing_pixels(obs; nside, nest = false, output_coordinate_system = Galactic,
                        detector_idx = :all, backend = CPUBackend()) -> PointingPixels

HEALPix pixel indices (1-based, as in Healpix.jl; subtract 1 for healpy),
orientation angles and HWP angles of every sample, ready for map-making. The
pixels are computed sample by sample, without storing `(θ, φ)`.
"""
function get_pointing_pixels(obs::Observation; nside::Integer, nest::Bool = false,
	output_coordinate_system::CoordinateSystem = Galactic, detector_idx = :all, backend::AbstractBackend = CPUBackend())
	plan, idx = _plan_data(obs, output_coordinate_system, detector_idx)
	resol = Resolution(nside)
	pixels, psi = if backend isa GPUBackend
		_require_gpu()
		_gpu_pointing_pixels(plan, obs.n_samples, length(idx), resol, nest, backend)
	else
		_cpu_pointing_pixels(plan, obs.n_samples, length(idx), resol, nest)
	end
	names = [obs.detectors[i].name for i in idx]
	return PointingPixels(pixels, psi, _hwp_angles(obs, plan), nside, nest, output_coordinate_system, names)
end

function _cpu_pointing_pixels(plan::_PlanData, nsamples::Int, ndet::Int, resol::Resolution, nest::Bool)
	pixels = Matrix{Int}(undef, nsamples, ndet)
	psi = Matrix{Float64}(undef, nsamples, ndet)
	bad = Threads.Atomic{Bool}(false)
	_foreach_block(plan, nsamples) do k_first, k_last, bp
		@threads for k in k_first:k_last
			@inbounds for n in 1:ndet
				θ, φ, ψ = compute_pointing_and_orientation(_sample_quat(bp, k, n))
				if !(isfinite(θ) & isfinite(φ) & isfinite(ψ))
					bad[] = true
					pixels[k, n] = 0
				else
					pixels[k, n] = _healpix_ang2pix(resol, nest, θ, φ)
				end
				psi[k, n] = ψ
			end
		end
	end
	bad[] && throw(ArgumentError(_NONFINITE_ANGLES))
	return pixels, psi
end

"""
    _foreach_chunk(f, plan, nsamples, ndet, resol, nest)

Call `f(first_sample, pixels, psi)` for consecutive blocks of samples, where
`pixels` and `psi` are `(block_length, ndet)` buffers.
"""
function _foreach_chunk(f::F, plan::_PlanData, nsamples::Int, ndet::Int, resol::Resolution, nest::Bool) where {F}
	len = _block_length(plan, nsamples)
	pixels = Matrix{Int}(undef, len, ndet)
	psi = Matrix{Float64}(undef, len, ndet)
	bad = Threads.Atomic{Bool}(false)
	_foreach_block(plan, nsamples) do first, last, bp
		m = last - first + 1
		@threads for j in 1:m
			@inbounds for n in 1:ndet
				θ, φ, ψ = compute_pointing_and_orientation(_sample_quat(bp, first + j - 1, n))
				if !(isfinite(θ) & isfinite(φ) & isfinite(ψ))
					bad[] = true
				else
					pixels[j, n] = _healpix_ang2pix(resol, nest, θ, φ)
					psi[j, n] = ψ
				end
			end
		end
		bad[] && throw(ArgumentError(_NONFINITE_ANGLES))
		f(first, view(pixels, 1:m, :), view(psi, 1:m, :))
	end
	return nothing
end

"""
    get_hitmap(obs; nside, output_coordinate_system = Galactic, nest = false,
               detector_idx = :all, backend = CPUBackend()) -> Vector{Int}

Number of samples falling in every HEALPix pixel, summed over the selected
detectors (and over the observations if `obs` is a vector), as `get_hitmap` in
litebird_sim.
"""
function get_hitmap(obs::Observation; kw...)
	return get_hitmap([obs]; kw...)
end

function get_hitmap(observations::AbstractVector{<:Observation}; nside::Integer,
	output_coordinate_system::CoordinateSystem = Galactic, nest::Bool = false, detector_idx = :all,
	backend::AbstractBackend = CPUBackend())
	backend isa GPUBackend && _require_gpu()
	resol = Resolution(nside)
	hits = zeros(Int, resol.numOfPixels)
	for obs in observations
		plan, idx = _plan_data(obs, output_coordinate_system, detector_idx)
		if backend isa GPUBackend
			_gpu_hits!(hits, plan, obs.n_samples, length(idx), resol, nest, backend)
		else
			_foreach_chunk(plan, obs.n_samples, length(idx), resol, nest) do _, pixels, _
				@inbounds for p in pixels
					hits[p] += 1
				end
			end
		end
	end
	return hits
end
