"""
    HMapsResult

Spin moments ("h-maps") of the scan, as returned by [`make_h_maps`](@ref):

    h_{n,m}(p) = (1 / N_p) Σ_{samples in p} exp(+i (n ψ + m χ))

with `ψ` the orientation angle (without the detector polarisation angle) and
`χ` the HWP angle, as `make_h_maps` in litebird_sim. For `(n, m) = (0, 0)` the
real part holds the hit count `N_p` instead; pixels without hits are 0.

Fields:
- `h`: `(n_couples, n_pixels, n_detectors)` ComplexF64 array (RING ordering).
- `n_m_couples`: the `(n, m)` pairs, in the order of the first axis of `h`.
- `detectors`: detector names, in the order of the last axis of `h`.
- `hits`: `(n_pixels, n_detectors)` hit counts, or `nothing` when read from a
  file without the `(0, 0)` couple.
- `nside`, `coordinate_system`, `duration_s` (the sum of the time spans of the
  observations; one value for all the detectors, as in litebird_sim),
  `sampling_rate_hz`, `detector_split` and `time_split` (always `"full"`).

Use [`h_map`](@ref) and [`hitmap`](@ref) to access a single map.
"""
struct HMapsResult
	h::Array{ComplexF64, 3}
	n_m_couples::Vector{Tuple{Int, Int}}
	detectors::Vector{String}
	hits::Union{Nothing, Matrix{Int}}
	nside::Int
	coordinate_system::CoordinateSystem
	duration_s::Float64
	sampling_rate_hz::Float64
	detector_split::String
	time_split::String
end

function Base.show(io::IO, res::HMapsResult)
	print(io, "HMapsResult(", length(res.detectors), " detectors, ", length(res.n_m_couples), " (n, m) couples, nside ",
		res.nside, ", ", res.coordinate_system, ")")
end

function _couple_index(res::HMapsResult, n::Integer, m::Integer)
	c = findfirst(==((Int(n), Int(m))), res.n_m_couples)
	c === nothing && throw(ArgumentError("(n, m) = ($n, $m) is not in the h-maps; available: $(res.n_m_couples)"))
	return c
end

function _detector_index(res::HMapsResult, det::AbstractString)
	d = findfirst(==(det), res.detectors)
	d === nothing && throw(ArgumentError("detector \"$det\" is not in the h-maps; available: $(res.detectors)"))
	return d
end

"""
    h_map(res::HMapsResult, det, n, m) -> AbstractVector{ComplexF64}

The map `h_{n,m}` of detector `det` (a view into `res.h`).
"""
h_map(res::HMapsResult, det::AbstractString, n::Integer, m::Integer) =
	view(res.h, _couple_index(res, n, m), :, _detector_index(res, det))

"""
    hitmap(res::HMapsResult, det) -> AbstractVector{Int}

Hit counts of detector `det`.
"""
function hitmap(res::HMapsResult, det::AbstractString)
	res.hits === nothing && throw(ArgumentError("these h-maps carry no hit counts (no (0, 0) couple)"))
	return view(res.hits, :, _detector_index(res, det))
end

"""
    nm_grid(spin_n, spin_m) -> Vector{Tuple{Int, Int}}

All couples `(n, m)` with `n ∈ spin_n` and `m ∈ spin_m`, e.g. the full grid
needed by SBM.
"""
nm_grid(spin_n, spin_m) = [(Int(n), Int(m)) for n in spin_n for m in spin_m]

function _normalize_couples(n_m_couples)
	couples = if n_m_couples isa AbstractMatrix
		size(n_m_couples, 2) == 2 || throw(ArgumentError("a matrix of (n, m) couples must have 2 columns"))
		[(n_m_couples[i, 1], n_m_couples[i, 2]) for i in axes(n_m_couples, 1)]
	else
		collect(n_m_couples)
	end
	isempty(couples) && throw(ArgumentError("no (n, m) couples given"))
	out = Tuple{Int, Int}[]
	for (n, m) in couples
		(isinteger(n) && isinteger(m)) || throw(ArgumentError("spin orders must be integers, got ($n, $m)"))
		push!(out, (Int(n), Int(m)))
	end
	allunique(out) || throw(ArgumentError("repeated (n, m) couples: $out"))
	return out
end

"""
    make_h_maps(observations; nside, n_m_couples = [(0, 0), (2, 0), (4, 0)],
                output_coordinate_system = Galactic, backend = CPUBackend()) -> HMapsResult

h-maps of every detector, as `make_h_maps` in litebird_sim (see
[`HMapsResult`](@ref)). `observations` is an [`Observation`](@ref) or a vector
of them; samples of the same detector (matched by name) are accumulated over
all observations before normalising. `n_m_couples` is a vector of integer
`(n, m)` tuples or an `(N, 2)` matrix; see [`nm_grid`](@ref).

Pointings are computed in chunks and never stored in full.
"""
make_h_maps(obs::Observation; kw...) = make_h_maps([obs]; kw...)

function make_h_maps(observations::AbstractVector{<:Observation}; nside::Integer,
	n_m_couples = [(0, 0), (2, 0), (4, 0)], output_coordinate_system::CoordinateSystem = Galactic,
	backend::AbstractBackend = CPUBackend())
	isempty(observations) && throw(ArgumentError("no observations given"))
	backend isa GPUBackend && _require_gpu()
	couples = _normalize_couples(n_m_couples)
	rates = unique(obs.sampling_rate_hz for obs in observations)
	length(rates) == 1 || throw(ArgumentError("all observations must have the same sampling rate, got $(rates) Hz"))

	names = String[]
	for obs in observations, det in obs.detectors
		det.name in names || push!(names, det.name)
	end
	resol = Resolution(nside)
	npix = resol.numOfPixels
	ncpl = length(couples)
	spin_n = Float64[n for (n, _) in couples]
	spin_m = Float64[m for (_, m) in couples]
	sums = zeros(ComplexF64, ncpl, npix, length(names))
	hits = zeros(Int, npix, length(names))

	for obs in observations
		plan, _ = _plan_data(obs, output_coordinate_system, :all)
		columns = Int[something(findfirst(==(det.name), names)) for det in obs.detectors]
		if backend isa GPUBackend
			_gpu_h_maps!(sums, hits, plan, obs.n_samples, columns, spin_n, spin_m, resol, backend)
			continue
		end
		χ = zeros(_block_length(plan, obs.n_samples))
		_foreach_chunk(plan, obs.n_samples, length(columns), resol, false) do first, pixels, psi
			for j in axes(pixels, 1)
				χ[j] = _hwp_at(plan, first + j - 1)
				isfinite(χ[j]) || throw(ArgumentError(_NONFINITE_ANGLES))
			end
			@inbounds for (n, col) in enumerate(columns), j in axes(pixels, 1)
				p = pixels[j, n]
				ψ = psi[j, n]
				x = χ[j]
				hits[p, col] += 1
				for c in 1:ncpl
					sums[c, p, col] += cis(spin_n[c] * ψ + spin_m[c] * x)
				end
			end
		end
	end

	c00 = findfirst(==((0, 0)), couples)
	@inbounds for col in axes(hits, 2), p in 1:npix
		nhit = hits[p, col]
		nhit == 0 && continue
		for c in 1:ncpl
			sums[c, p, col] = c == c00 ? complex(Float64(nhit)) : sums[c, p, col] / nhit
		end
	end
	duration_s = sum(get_time_span, observations)
	return HMapsResult(sums, couples, names, hits, nside, output_coordinate_system, duration_s, rates[1], "full", "full")
end

"""
    save_h_maps(res::HMapsResult, output_directory) -> Vector{String}

Write one HDF5 file per detector, `h_maps_det_<name>.h5`, in the layout used by
`save_h_maps` in litebird_sim: root attributes `coordinate_system`, `det`,
`detector_split`, `time_split`, `duration_s` and `sampling_rate_Hz`, and one
group `"n,m"` per couple holding the datasets `Re` and `Im` (RING ordering).
The `(0, 0)` group holds the hit counts, and is written also when `(0, 0)` is
not among the couples (SBM and the map-makers need the hits). Returns the paths
of the files written.

The string attributes are scalar variable-length UTF-8 strings, as written by
h5py, so that litebird_sim reads them as `str`.
"""
function save_h_maps(res::HMapsResult, output_directory::AbstractString)
	mkpath(output_directory)
	paths = String[]
	for (d, det) in enumerate(res.detectors)
		path = joinpath(output_directory, "h_maps_det_$(det).h5")
		h5open(path, "w") do f
			_write_vlen_string_attribute(f, "coordinate_system", string(res.coordinate_system))
			_write_vlen_string_attribute(f, "det", det)
			_write_vlen_string_attribute(f, "detector_split", res.detector_split)
			_write_vlen_string_attribute(f, "time_split", res.time_split)
			attrs = attributes(f)
			attrs["duration_s"] = res.duration_s
			attrs["sampling_rate_Hz"] = res.sampling_rate_hz
			for (c, (n, m)) in enumerate(res.n_m_couples)
				f["$n,$m/Re"] = real.(res.h[c, :, d])
				f["$n,$m/Im"] = imag.(res.h[c, :, d])
			end
			if (0, 0) ∉ res.n_m_couples && res.hits !== nothing
				f["0,0/Re"] = Float64.(res.hits[:, d])
				f["0,0/Im"] = zeros(size(res.hits, 1))
			end
		end
		push!(paths, path)
	end
	return paths
end

# HDF5.jl writes a `String` attribute as a fixed-length string, which h5py
# returns as `bytes`; litebird_sim needs `str` (it looks up the coordinate
# system by name).
function _write_vlen_string_attribute(obj, name::AbstractString, value::AbstractString)
	dtype = HDF5.datatype(String)  # variable length, UTF-8
	dspace = HDF5.Dataspace(HDF5.API.h5s_create(HDF5.API.H5S_SCALAR))
	attr = create_attribute(obj, name, dtype, dspace)
	try
		write_attribute(attr, dtype, String(value))
	finally
		close(attr)
		close(dspace)
		close(dtype)
	end
	return nothing
end

"""
    load_h_maps_from_file(path) -> HMapsResult

Read the h-maps of one detector written by [`save_h_maps`](@ref) or by
`save_h_maps` in litebird_sim.
"""
function load_h_maps_from_file(path::AbstractString)
	h5open(path, "r") do f
		attrs = attributes(f)
		coord_name = read(attrs["coordinate_system"])
		coord = coord_name == "Galactic" ? Galactic : coord_name == "Ecliptic" ? Ecliptic :
				throw(ArgumentError("unknown coordinate system \"$coord_name\" in $path"))
		couples = Tuple{Int, Int}[]
		maps = Vector{Vector{ComplexF64}}()
		for key in keys(f)
			parts = split(key, ",")
			length(parts) == 2 || continue
			push!(couples, (parse(Int, parts[1]), parse(Int, parts[2])))
			push!(maps, complex.(read(f[key]["Re"]), read(f[key]["Im"])))
		end
		isempty(couples) && throw(ArgumentError("no h-maps found in $path"))
		npix = length(maps[1])
		h = Array{ComplexF64}(undef, length(couples), npix, 1)
		for c in eachindex(maps)
			h[c, :, 1] = maps[c]
		end
		c00 = findfirst(==((0, 0)), couples)
		hits = c00 === nothing ? nothing : reshape(round.(Int, real.(maps[c00])), npix, 1)
		return HMapsResult(h, couples, [String(read(attrs["det"]))], hits, npix2nside(npix), coord,
			Float64(read(attrs["duration_s"])), Float64(read(attrs["sampling_rate_Hz"])),
			String(read(attrs["detector_split"])), String(read(attrs["time_split"])))
	end
end
