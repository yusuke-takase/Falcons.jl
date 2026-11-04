"""
    Observation(attitude, instrument, detectors; start_time_s, duration_s, hwp = nothing,
                start_time_simulation_s = start_time_s, spin_perturbation = nothing,
                delta_time_s = nothing)

A continuous stretch of observation by a set of detectors, the counterpart of
an `Observation` in litebird_sim without the TOD.

- `attitude`: the rotation from the spin-axis frame to Ecliptic coordinates.
  Either a [`SpinningScanningStrategy`](@ref) (evaluated analytically) or a
  time-dependent [`RotQuaternion`](@ref), e.g. a measured attitude converted to
  this frame by the user.
- `instrument`: an [`InstrumentInfo`](@ref) giving the boresight-to-spin rotation.
- `detectors`: a [`DetectorInfo`](@ref) or a vector of them. All detectors must
  share the same `sampling_rate_hz`.
- `start_time_s`: absolute time of the first sample (s). Spin, precession and
  the anti-Sun direction are evaluated at absolute times, as litebird_sim does
  for a floating-point start time.
- `duration_s`: the observation has `floor(sampling_rate_hz * duration_s)`
  samples at `start_time_s + k / sampling_rate_hz`.
- `hwp`: an [`IdealHWP`](@ref), or `nothing`.
- `start_time_simulation_s`: origin of the HWP phase. Observations split from
  the same simulation share it, so the HWP angle is continuous across them.
- `spin_perturbation`: optional [`RotQuaternion`](@ref) multiplied on the right
  of the spin-to-Ecliptic rotation, i.e. a rotation expressed in the spin-axis
  frame (the position of `SpacecraftCoord` in litebird_sim).
- `delta_time_s`: `nothing` computes every sample exactly. A number reproduces
  litebird_sim: the boresight rotation is computed on anchors every
  `delta_time_s` seconds from `start_time_simulation_s` and interpolated with
  slerp.

The detector quaternion at a sample is

    q(t) = A(t) · S(t) · bore2spin · d(t)

with `A` the attitude, `S` the spin perturbation and `d` the detector
quaternion; time series are interpolated to `t`.
"""
struct Observation{A <: Union{SpinningScanningStrategy, RotQuaternion}, H <: Union{Nothing, IdealHWP}, P <: Union{Nothing, RotQuaternion}}
	attitude::A
	instrument::InstrumentInfo
	detectors::Vector{DetectorInfo}
	start_time_s::Float64
	n_samples::Int
	sampling_rate_hz::Float64
	hwp::H
	start_time_simulation_s::Float64
	spin_perturbation::P
	delta_time_s::Union{Nothing, Float64}

	function Observation(attitude::A, instrument::InstrumentInfo, detectors::Vector{DetectorInfo}, start_time_s::Real,
		n_samples::Integer, hwp::H, start_time_simulation_s::Real, spin_perturbation::P, delta_time_s) where {A, H, P}
		isempty(detectors) && throw(ArgumentError("an Observation needs at least one detector"))
		rates = unique(d.sampling_rate_hz for d in detectors)
		length(rates) == 1 || throw(ArgumentError("all detectors must have the same sampling rate, got $(rates) Hz"))
		rates[1] > 0 || throw(ArgumentError("the detectors' sampling_rate_hz must be positive"))
		n_samples >= 1 || throw(ArgumentError("the observation has no samples"))
		if delta_time_s !== nothing && !(delta_time_s > 0)
			throw(ArgumentError("delta_time_s must be positive or nothing"))
		end
		names = [d.name for d in detectors]
		allunique(names) || throw(ArgumentError("detector names must be unique, got $(names)"))
		return new{A, H, P}(attitude, instrument, detectors, start_time_s, n_samples, rates[1], hwp,
			start_time_simulation_s, spin_perturbation, delta_time_s === nothing ? nothing : Float64(delta_time_s))
	end
end

function Observation(attitude, instrument::InstrumentInfo, detectors; start_time_s::Real, duration_s::Real,
	hwp = nothing, start_time_simulation_s::Real = start_time_s, spin_perturbation = nothing, delta_time_s = nothing)
	dets = detectors isa DetectorInfo ? [detectors] : collect(DetectorInfo, detectors)
	isempty(dets) && throw(ArgumentError("an Observation needs at least one detector"))
	n_samples = floor(Int, dets[1].sampling_rate_hz * duration_s)
	return Observation(attitude, instrument, dets, start_time_s, n_samples, hwp, start_time_simulation_s,
		spin_perturbation, delta_time_s)
end

"Number of detectors of `obs`."
n_detectors(obs::Observation) = length(obs.detectors)

"Time span covered by the samples of `obs` (s), `n_samples / sampling_rate_hz`."
get_time_span(obs::Observation) = obs.n_samples / obs.sampling_rate_hz

"""
    sample_time(obs, k) -> t

Absolute time (s) of the `k`-th sample (`k = 1, …, n_samples`).
"""
@inline sample_time(obs::Observation, k::Integer) = obs.start_time_s + (k - 1) * (1 / obs.sampling_rate_hz)

"""
    get_times(obs) -> Vector{Float64}

Absolute times (s) of all samples of `obs`.
"""
get_times(obs::Observation) = [sample_time(obs, k) for k in 1:obs.n_samples]

"""
    split_observation(obs, n) -> Vector{Observation}

Split `obs` into `n` consecutive observations, distributing the samples as
evenly as possible (the first observations get one more sample when the
number of samples is not a multiple of `n`), as litebird_sim does when it
creates several observations per detector. The HWP origin
`start_time_simulation_s` is kept, so the HWP angle is continuous.
"""
function split_observation(obs::Observation, n::Integer)
	n >= 1 || throw(ArgumentError("n must be ≥ 1, got $n"))
	n <= obs.n_samples || throw(ArgumentError("cannot split $(obs.n_samples) samples into $n observations"))
	base, leftovers = divrem(obs.n_samples, n)
	result = typeof(obs)[]
	cur_time = obs.start_time_s
	for i in 1:n
		nsamples = base + (i <= leftovers ? 1 : 0)
		push!(result, Observation(obs.attitude, obs.instrument, obs.detectors, cur_time, nsamples, obs.hwp,
			obs.start_time_simulation_s, obs.spin_perturbation, obs.delta_time_s))
		cur_time += nsamples / obs.sampling_rate_hz
	end
	return result
end

function Base.show(io::IO, obs::Observation)
	print(io, "Observation(", n_detectors(obs), " detectors, ", obs.n_samples, " samples at ", obs.sampling_rate_hz,
		" Hz from t = ", obs.start_time_s, " s)")
end
