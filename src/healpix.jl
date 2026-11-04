# HEALPix pixelisation that can run inside GPU kernels: the same operations as
# `ang2pixRing` / `ang2pixNest` in Healpix.jl, without exceptions, table
# lookups or `mod2pi`. Both expect θ ∈ [0, π] and φ ∈ [0, 2π) (as returned by
# `compute_pointing_and_orientation`) and return 1-based pixel indices.

@inline _cachesintheta(θ) = ((θ < 0.01) || (θ > 3.14159 - 0.01)) ? (sin(θ), true) : (0.0, false)

@inline _ifloor(x) = unsafe_trunc(Int, floor(x))

"Interleave the bits of `v` with zeros (Healpix.jl's `spreadbits`)."
@inline function _spread_bits(v::Int)
	x = UInt64(v) & 0x00000000ffffffff
	x = (x | (x << 16)) & 0x0000ffff0000ffff
	x = (x | (x << 8)) & 0x00ff00ff00ff00ff
	x = (x | (x << 4)) & 0x0f0f0f0f0f0f0f0f
	x = (x | (x << 2)) & 0x3333333333333333
	x = (x | (x << 1)) & 0x5555555555555555
	return Int(x)
end

@inline _xyf2pix_nest(r::Resolution, ix, iy, face) = face << (2 * r.order) + _spread_bits(ix) + 2 * _spread_bits(iy) + 1

@inline function _ang2pix_ring(r::Resolution, θ, φ)
	z = cos(θ)
	z_abs = abs(z)
	tt = mod(φ * 2 / π, 4.0e0)
	sintheta, havesintheta = _cachesintheta(θ)
	if 3 * z_abs <= 2
		nl4 = 4 * r.nside
		temp1 = r.nside * (0.5 + tt)
		temp2 = r.nside * z * 0.75
		jp = _ifloor(temp1 - temp2)
		jm = _ifloor(temp1 + temp2)
		ir = r.nside + 1 + jp - jm
		kshift = 1 - (ir & 1)
		t1 = jp + jm - r.nside + kshift + 1 + 2nl4
		ip = (r.order > 0) ? ((t1 >> 1) & (nl4 - 1)) : mod(t1 >> 1, nl4)
		return r.ncap + (ir - 1) * nl4 + ip + 1
	end
	tp = tt - _ifloor(tt)
	tmp = (z_abs < 0.99 || !havesintheta) ? r.nside * sqrt(3 * (1 - z_abs)) : r.nside * sintheta / sqrt((1 + z_abs) / 3)
	jp = _ifloor(tp * tmp)
	jm = _ifloor((1 - tp) * tmp)
	ir = jp + jm + 1
	ip = _ifloor(tt * ir)
	return z > 0 ? 2 * ir * (ir - 1) + ip + 1 : r.numOfPixels - 2 * ir * (ir + 1) + ip + 1
end

@inline function _ang2pix_nest(r::Resolution, θ, φ)
	z = cos(θ)
	z_abs = abs(z)
	scaled_phi = φ / (π / 2)
	sintheta, havesintheta = _cachesintheta(θ)
	if 3 * z_abs <= 2
		temp1 = r.nside * (0.5 + scaled_phi)
		temp2 = r.nside * (z * 0.75)
		jp = _ifloor(temp1 - temp2)
		jm = _ifloor(temp1 + temp2)
		ifp = jp >> r.order
		ifm = jm >> r.order
		face = ifp == ifm ? (ifp | 4) : ifp < ifm ? ifp : ifm + 8
		ix = mod(jm, r.nside)
		iy = r.nside - mod(jp, r.nside) - 1
		return _xyf2pix_nest(r, ix, iy, face)
	end
	ntt = _ifloor(scaled_phi)
	ntt >= 4 && (ntt = 3)
	tp = scaled_phi - ntt
	tmp = (z_abs < 0.99 || !havesintheta) ? r.nside * sqrt(3 * (1 - z_abs)) : r.nside * sintheta / sqrt((1 + z_abs) / 3)
	jp = min(_ifloor(tp * tmp), r.nside - 1)
	jm = min(_ifloor((1 - tp) * tmp), r.nside - 1)
	return z >= 0 ? _xyf2pix_nest(r, r.nside - jm - 1, r.nside - jp - 1, ntt) : _xyf2pix_nest(r, jp, jm, ntt + 8)
end

@inline _ang2pix(r::Resolution, nest::Bool, θ, φ) = nest ? _ang2pix_nest(r, θ, φ) : _ang2pix_ring(r, θ, φ)

# 2π split for an accurate reduction of large angles (Cody–Waite with fma).
const _TWO_PI_HI = 6.283185307179586
const _TWO_PI_LO = 2.4492935982947064e-16

"`x` reduced to [−π, π], accurate for |x| up to ~1e15 and usable in GPU kernels."
@inline function _rem2pi(x::Float64)
	k = round(x / _TWO_PI_HI)
	r = fma(-k, _TWO_PI_HI, x)
	return fma(-k, _TWO_PI_LO, r)
end
