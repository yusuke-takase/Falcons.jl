"""
CUDA backend of Falcons, loaded by `using CUDA`.

The kernels evaluate the same `_sample_quat` as the CPU code on a copy of the
`_PlanData` whose arrays live on the device, so both backends follow the same
operations in Float64. Samples are split into contiguous ranges, one per GPU;
hit maps and h-maps are accumulated with atomics and summed on the host. In
interpolated mode the anchors of each block of samples are computed on the
host and uploaded before the kernels of that block.
"""
module FalconsCUDAExt

using CUDA
using Falcons
using Falcons: _PlanData, _Series, _DetectorQuats, _sample_quat, _hwp_at, _ang2pix, _ang2pix_ring,
	_anchor_block, _block_length, _with_anchors, _NONFINITE_ANGLES, compute_pointing_and_orientation, GPUBackend,
	Quaternion
using Healpix: Resolution

const THREADS = 256
const MAX_BLOCKS = 65535

_blocks(n) = max(1, min(cld(n, THREADS), MAX_BLOCKS))

function _device_ids(backend::GPUBackend)
	CUDA.functional() || throw(ArgumentError("GPUBackend: CUDA is not functional on this machine"))
	available = length(CUDA.devices())
	ids = backend.gpu_ids === nothing ? collect(0:backend.num_gpus-1) : backend.gpu_ids
	for id in ids
		0 <= id < available || throw(ArgumentError("GPU id $id is not available (this machine has $available GPUs)"))
	end
	return ids
end

"Split `1:n` into `parts` contiguous ranges (the first ones one longer)."
function _ranges(n::Int, parts::Int)
	base, extra = divrem(n, parts)
	out = UnitRange{Int}[]
	first = 1
	for i in 1:parts
		len = base + (i <= extra ? 1 : 0)
		len > 0 && push!(out, first:first+len-1)
		first += len
	end
	return out
end

"Run `f(range)` on every selected GPU, each on its share of the samples."
function _on_gpus(f, backend::GPUBackend, nsamples::Int)
	ids = _device_ids(backend)
	ranges = _ranges(nsamples, length(ids))
	tasks = map(zip(ids, ranges)) do (id, range)
		Threads.@spawn begin
			CUDA.device!(id)
			result = f(range)
			CUDA.synchronize()
			result
		end
	end
	return fetch.(tasks)
end

# Copy the arrays of a plan to the current device. The CuArrays are returned
# as well so that they stay alive while the kernels run.
function _upload(p::_PlanData)
	keep = CuArray[]
	up(x) = (c = CuArray(x); push!(keep, c); cudaconvert(c))
	series(s::_Series) = _Series(up(s.quats), s.t0, s.rate)
	series(::Nothing) = nothing
	attitude = p.attitude isa _Series ? series(p.attitude) : p.attitude
	d = p.dets
	dets = _DetectorQuats(up(d.constant), up(d.series), up(d.offset), up(d.len), up(d.t0), up(d.rate))
	anchors = up(fill(one(Quaternion{Float64}), 2, 1))  # placeholder: see _device_block
	dp = _PlanData(attitude, series(p.perturbation), p.b, p.g, p.galactic, dets, p.interp, anchors, 0, p.n_anchors,
		p.anchor_j0, p.anchor_tsim, p.anchor_dt, p.anchor_t0, p.anchor_rate, p.start_time_s, p.dt, p.has_hwp,
		p.hwp_omega, p.hwp_start, p.hwp_dt0)
	return dp, keep
end

"""
    _device_block(dp, plan, k_first, k_last) -> (device plan, anchors or nothing)

The device plan `dp` with the anchors of the samples `k_first:k_last`,
computed on the host from `plan` and uploaded. The returned CuArray must be
kept alive while the kernels of the block run (and can be freed afterwards).
In exact mode `dp` is returned unchanged.
"""
function _device_block(dp::_PlanData, plan::_PlanData, k_first::Int, k_last::Int)
	plan.interp || return dp, nothing
	block = _anchor_block(plan, k_first, k_last)
	anchors = CuArray(block.anchors)
	return _with_anchors(dp, cudaconvert(anchors), block.anchor_first), anchors
end

_free!(::Nothing) = nothing
_free!(a::CuArray) = CUDA.unsafe_free!(a)

@inline function _index(m)
	i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
	return i, gridDim().x * blockDim().x
end

function _k_pointings!(out, p, first, m, ndet)
	i, stride = _index(m)
	@inbounds while i <= m * ndet
		j = (i - 1) % m + 1
		n = (i - 1) ÷ m + 1
		θ, φ, ψ = compute_pointing_and_orientation(_sample_quat(p, first + j - 1, n))
		out[1, j, n] = θ
		out[2, j, n] = φ
		out[3, j, n] = ψ
		i += stride
	end
	return nothing
end

# The kernels that index maps skip samples with non-finite angles (whose pixel
# index would be undefined) and set `bad[1]`; the host then throws, as the CPU
# backend does. The inputs are checked on the host first, so this is a backstop.

function _k_pixels!(pixels, psi, bad, p, resol, nest, first, m, ndet)
	i, stride = _index(m)
	@inbounds while i <= m * ndet
		j = (i - 1) % m + 1
		n = (i - 1) ÷ m + 1
		θ, φ, ψ = compute_pointing_and_orientation(_sample_quat(p, first + j - 1, n))
		if !(isfinite(θ) & isfinite(φ) & isfinite(ψ))
			bad[1] = Int32(1)
			pixels[j, n] = 0
		else
			pixels[j, n] = _ang2pix(resol, nest, θ, φ)
		end
		psi[j, n] = ψ
		i += stride
	end
	return nothing
end

function _k_hits!(hits, bad, p, resol, nest, first, m, ndet)
	i, stride = _index(m)
	@inbounds while i <= m * ndet
		j = (i - 1) % m + 1
		n = (i - 1) ÷ m + 1
		θ, φ, ψ = compute_pointing_and_orientation(_sample_quat(p, first + j - 1, n))
		if !(isfinite(θ) & isfinite(φ) & isfinite(ψ))
			bad[1] = Int32(1)
		else
			CUDA.@atomic hits[_ang2pix(resol, nest, θ, φ)] += 1
		end
		i += stride
	end
	return nothing
end

function _k_h_maps!(re, im, hits, bad, p, resol, first, m, ndet, columns, spin_n, spin_m)
	i, stride = _index(m)
	ncpl = length(spin_n)
	@inbounds while i <= m * ndet
		j = (i - 1) % m + 1
		n = (i - 1) ÷ m + 1
		k = first + j - 1
		θ, φ, ψ = compute_pointing_and_orientation(_sample_quat(p, k, n))
		χ = _hwp_at(p, k)
		if !(isfinite(θ) & isfinite(φ) & isfinite(ψ) & isfinite(χ))
			bad[1] = Int32(1)
			i += stride
			continue
		end
		pix = _ang2pix_ring(resol, θ, φ)
		col = columns[n]
		CUDA.@atomic hits[pix, col] += 1
		for c in 1:ncpl
			s, co = sincos(spin_n[c] * ψ + spin_m[c] * χ)
			CUDA.@atomic re[c, pix, col] += co
			CUDA.@atomic im[c, pix, col] += s
		end
		i += stride
	end
	return nothing
end

# Samples per launch for the kernels that write per-sample outputs.
_chunk(ndet) = max(1, (1 << 24) ÷ ndet)

# Samples per kernel launch in interpolated mode (one block of anchors each).
_block(plan, range, ndet) = _block_length(plan, length(range), _chunk(ndet))

function Falcons._gpu_pointings(plan::_PlanData, nsamples::Int, ndet::Int, backend::GPUBackend)
	pt = Array{Float64}(undef, 3, nsamples, ndet)
	_on_gpus(backend, nsamples) do range
		dp, _keep = _upload(plan)
		len = _block(plan, range, ndet)
		for first in range.start:len:range.stop
			m = min(len, range.stop - first + 1)
			bp, anchors = _device_block(dp, plan, first, first + m - 1)
			out = CuArray{Float64}(undef, 3, m, ndet)
			@cuda threads = THREADS blocks = _blocks(m * ndet) _k_pointings!(out, bp, first, m, ndet)
			pt[:, first:first+m-1, :] = Array(out)
			CUDA.unsafe_free!(out)
			_free!(anchors)
		end
		nothing
	end
	return pt
end

function Falcons._gpu_pointing_pixels(plan::_PlanData, nsamples::Int, ndet::Int, resol::Resolution, nest::Bool,
	backend::GPUBackend)
	pixels = Matrix{Int}(undef, nsamples, ndet)
	psi = Matrix{Float64}(undef, nsamples, ndet)
	bad = _on_gpus(backend, nsamples) do range
		dp, _keep = _upload(plan)
		d_bad = CUDA.zeros(Int32, 1)
		len = _block(plan, range, ndet)
		for first in range.start:len:range.stop
			m = min(len, range.stop - first + 1)
			bp, anchors = _device_block(dp, plan, first, first + m - 1)
			d_pix = CuArray{Int}(undef, m, ndet)
			d_psi = CuArray{Float64}(undef, m, ndet)
			@cuda threads = THREADS blocks = _blocks(m * ndet) _k_pixels!(d_pix, d_psi, d_bad, bp, resol, nest, first, m,
				ndet)
			pixels[first:first+m-1, :] = Array(d_pix)
			psi[first:first+m-1, :] = Array(d_psi)
			CUDA.unsafe_free!(d_pix)
			CUDA.unsafe_free!(d_psi)
			_free!(anchors)
		end
		Array(d_bad)[1] != 0
	end
	_check_finite_angles(bad)
	return pixels, psi
end

_check_finite_angles(bad) = any(bad) ? throw(ArgumentError(_NONFINITE_ANGLES)) : nothing

function Falcons._gpu_hits!(hits::Vector{Int}, plan::_PlanData, nsamples::Int, ndet::Int, resol::Resolution, nest::Bool,
	backend::GPUBackend)
	partials = _on_gpus(backend, nsamples) do range
		dp, _keep = _upload(plan)
		d_hits = CUDA.zeros(Int, length(hits))
		d_bad = CUDA.zeros(Int32, 1)
		len = plan.interp ? _block(plan, range, ndet) : length(range)
		for first in range.start:len:range.stop
			m = min(len, range.stop - first + 1)
			bp, anchors = _device_block(dp, plan, first, first + m - 1)
			@cuda threads = THREADS blocks = _blocks(m * ndet) _k_hits!(d_hits, d_bad, bp, resol, nest, first, m, ndet)
			_free!(anchors)
		end
		(Array(d_hits), Array(d_bad)[1] != 0)
	end
	_check_finite_angles(last.(partials))
	for (part, _) in partials
		hits .+= part
	end
	return hits
end

# Function barrier: the partial results come out of `fetch` untyped.
function _add_complex!(sums::Array{ComplexF64, 3}, re::Array{Float64, 3}, im::Array{Float64, 3})
	@inbounds for i in eachindex(sums, re, im)
		sums[i] += complex(re[i], im[i])
	end
	return sums
end

function Falcons._gpu_h_maps!(sums::Array{ComplexF64, 3}, hits::Matrix{Int}, plan::_PlanData, nsamples::Int,
	columns::Vector{Int}, spin_n::Vector{Float64}, spin_m::Vector{Float64}, resol::Resolution, backend::GPUBackend)
	ndet = length(columns)
	ids = _device_ids(backend)
	partials = _on_gpus(backend, nsamples) do range
		dp, _keep = _upload(plan)
		d_re = CUDA.zeros(Float64, size(sums))
		d_im = CUDA.zeros(Float64, size(sums))
		d_hits = CUDA.zeros(Int, size(hits))
		d_bad = CUDA.zeros(Int32, 1)
		d_cols, d_n, d_m = CuArray(columns), CuArray(spin_n), CuArray(spin_m)
		len = plan.interp ? _block(plan, range, ndet) : length(range)
		for first in range.start:len:range.stop
			m = min(len, range.stop - first + 1)
			bp, anchors = _device_block(dp, plan, first, first + m - 1)
			@cuda threads = THREADS blocks = _blocks(m * ndet) _k_h_maps!(d_re, d_im, d_hits, d_bad, bp, resol, first, m,
				ndet, d_cols, d_n, d_m)
			_free!(anchors)
		end
		(d_re, d_im, d_hits, Array(d_bad)[1] != 0)
	end
	# Reduce the partial sums on the first device (device-to-device copies),
	# then download once. `device!` with a function restores the caller's
	# device afterwards.
	CUDA.device!(ids[1]) do
		try
			_check_finite_angles(last.(partials))
			re, im, h = partials[1]
			for (o_re, o_im, o_h) in partials[2:end]
				re .+= CuArray(o_re)
				im .+= CuArray(o_im)
				h .+= CuArray(o_h)
			end
			_add_complex!(sums, Array(re), Array(im))
			hits .+= Array(h)
		finally
			for part in partials, a in part[1:3]
				CUDA.unsafe_free!(a)
			end
		end
	end
	return sums, hits
end

end
