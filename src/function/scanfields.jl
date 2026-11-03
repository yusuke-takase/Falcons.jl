"""
    ScanField

Hit map and spin moments of one scan, as returned by [`get_scanfield`](@ref).

- `hitmap`: number of samples per RING pixel.
- `h`: `h[i, j, p] = ⟨exp(-i(n ψ + m ϕ))⟩` in pixel `p` for `n = n[i]`,
  `m = m[j]`, with `ψ` the orientation angle and `ϕ` the HWP angle. Pixels
  without hits are `NaN`.
- `quantify`: summary statistics of `|h|²` (see [`get_hnm_quantify`](@ref)).
- `n`, `m`: the spin orders.
- `ss`: the scanning strategy that produced the field.
"""
struct ScanField{S <: ScanningStrategy}
	hitmap::Vector{Int64}
	h::Array{ComplexF64, 3}
	quantify::DataFrame
	n::Vector
	m::Vector
	ss::S
end

function h_nm(field::ScanField, n::Real, m::Real)
	return h_nm(field.h, field.n, field.m, n, m)
end

function h_nm(hₙₘ::Array, spin_n::AbstractVector, spin_m::AbstractVector, n::Real, m::Real)
	i = findfirst(==(n), spin_n)
	j = findfirst(==(m), spin_m)
	(i === nothing || j === nothing) && throw(ArgumentError("(n, m) = ($n, $m) is not in the scan field"))
	return hₙₘ[i, j, :]
end

_nanmean(x) = mean(filter(!isnan, x))
_nanstd(x) = std(filter(!isnan, x))

function get_hnm_quantify(hₙₘ, spin_n, spin_m)
	df = DataFrame(n            = repeat(spin_n, inner=length(spin_m)),
	               m            = repeat(spin_m, outer=length(spin_n)),
	               mean         = zeros(length(spin_n)*length(spin_m)),
	               std          = zeros(length(spin_n)*length(spin_m)),
	               nanmean      = zeros(length(spin_n)*length(spin_m)),
	               nanstd       = zeros(length(spin_n)*length(spin_m)),
	               nan2one_mean = zeros(length(spin_n)*length(spin_m)),
	               nan2one_std  = zeros(length(spin_n)*length(spin_m)),
	)
	for _n in spin_n
		for _m in spin_m
			row = findfirst((r -> r.n == _n && r.m == _m), eachrow(df))
			abs2_field = abs2.(h_nm(hₙₘ, spin_n, spin_m, _n, _m))
			df[row, :mean]         = mean(abs2_field)
			df[row, :std]          = std(abs2_field)
			df[row, :nanmean]      = _nanmean(abs2_field)
			df[row, :nanstd]       = _nanstd(abs2_field)
			df[row, :nan2one_mean] = mean(replace(abs2_field, NaN=>1.0))
			df[row, :nan2one_std]  = std(replace(abs2_field, NaN=>1.0))
		end
	end
	return df
end

"""
    _scanfield(ss, pointings_between; division, spin_n, spin_m, verbose = false)

Accumulate the hit map and spin moments, computing the pointings of
`division` consecutive time chunks with `pointings_between(start, stop)`.
Chunk boundaries are whole seconds, so `division` need not divide the
duration. With `verbose = true` a message is logged after every chunk.
"""
function _scanfield(ss::ScanningStrategy, pointings_between::F; division::Int, spin_n::AbstractVector, spin_m::AbstractVector, verbose::Bool = false) where {F}
	division >= 1 || throw(ArgumentError("division must be ≥ 1, got $division"))
	resol  = Resolution(ss.nside)
	npix   = nside2npix(ss.nside)
	ω_hwp  = rpm2angfreq(ss.hwp_rpm)
	hitmap = zeros(Int64, npix)
	hₙₘ    = zeros(ComplexF64, (length(spin_n), length(spin_m), npix))
	polang = [get_pol_angle(ss, j) for j in eachindex(ss.quat)]
	for i in 1:division
		chunk_start = div((i - 1) * ss.duration, division)
		chunk_stop  = div(i * ss.duration, division)
		theta, phi, psi, time = pointings_between(chunk_start, chunk_stop)
		@inbounds for j in axes(theta, 2), k in eachindex(time)
			Ω = ang2pixRing(resol, theta[k, j], phi[k, j])
			ψ = psi[k, j]
			ϕ = mod2pi(ω_hwp * time[k]) + polang[j]
			hitmap[Ω] += 1
			for _m in eachindex(spin_m), _n in eachindex(spin_n)
				hₙₘ[_n, _m, Ω] += cis(-(spin_n[_n] * ψ + spin_m[_m] * ϕ))
			end
		end
		verbose && @info "scan field: chunk $i/$division done" chunk_stop
	end
	for p in 1:npix, _m in eachindex(spin_m), _n in eachindex(spin_n)
		hₙₘ[_n, _m, p] /= hitmap[p]
	end
	df = get_hnm_quantify(hₙₘ, spin_n, spin_m)
	return ScanField(hitmap, hₙₘ, df, collect(spin_n), collect(spin_m), ss)
end

"""
    get_scanfield(ss; division, spin_n, spin_m, verbose = false) -> ScanField

Hit map and spin moments `h_{n,m}` of the whole scan (`0 ≤ t < ss.duration`),
for every combination of `n ∈ spin_n` and `m ∈ spin_m`. The pointings are
computed in `division` time chunks to bound memory use.
"""
function get_scanfield(ss::ScanningStrategy; division::Int, spin_n::AbstractVector, spin_m::AbstractVector, verbose::Bool = false)
	return _scanfield(ss, (a, b) -> get_pointings(ss, a, b); division, spin_n, spin_m, verbose)
end

"""
    get_scanfield_step(ss; division, spin_n, spin_m, step, verbose = false) -> ScanField

Same as [`get_scanfield`](@ref), with the anti-Sun direction advancing every
`step` seconds (see `get_pointings(ss, start, stop, step)`).
"""
function get_scanfield_step(ss::ScanningStrategy; division::Int, spin_n::AbstractVector, spin_m::AbstractVector, step, verbose::Bool = false)
	return _scanfield(ss, (a, b) -> get_pointings(ss, a, b, step); division, spin_n, spin_m, verbose)
end
