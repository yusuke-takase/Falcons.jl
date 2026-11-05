# GPU backend

[`get_pointings`](@ref), [`get_pointing_pixels`](@ref), [`get_hitmap`](@ref) and [`make_h_maps`](@ref)
accept a `backend` keyword. Pass a [`GPUBackend`](@ref) to compute on CUDA GPUs:

```julia
using CUDA        # loads the FalconsCUDAExt extension
using Falcons

res = make_h_maps(obs; nside = 512, n_m_couples = nm_grid(0:4, (-4, 0, 4)), backend = GPUBackend())
hits = get_hitmap(obs; nside = 512, backend = GPUBackend(num_gpus = 4))
```

- `GPUBackend(num_gpus = N)` uses the first `N` devices and `GPUBackend(gpu_ids = [...])` selects devices by
  their 0-based ids. The samples are split between the devices and the partial maps are summed. Start Julia
  with at least as many threads as GPUs (`julia -t 4`) so that the devices run concurrently.
- All modes are supported: analytic or time-series attitude, spin perturbations, time-dependent detector
  quaternions, exact or interpolated pointings, Ecliptic or Galactic coordinates, RING or NEST pixels.
- Computations are in Float64 and follow the same operations as the CPU backend, so pointings agree to
  ~1e-12 rad, pixel indices and hit maps exactly, and h-maps to ~1e-10 (accumulation order).
- The inputs are checked on the host before any kernel runs: a zero or non-finite quaternion (for instance
  a missing sample of an attitude time series) or a non-finite time or angle raises an `ArgumentError`, as
  on the CPU.

CUDA is a weak dependency: Falcons works without it, and calling a GPU backend without `using CUDA` raises
an `ArgumentError`.

## Memory

Each GPU accumulates the full maps of all detectors: [`make_h_maps`](@ref) allocates, on every device, the
real and imaginary parts of the sums, `2 × n_couples × n_pixels × n_detectors` Float64 values, and the hit
counts, `n_pixels × n_detectors` integers, that is about `(16 n_couples + 8) × n_pixels` bytes per detector.
Splitting the samples between several GPUs makes the computation faster but does not reduce this memory:
at nside 512 with 15 couples a detector takes about 0.78 GB, so about 20 detectors fill a 16 GB device
whatever the number of GPUs. The host holds the same sums in ComplexF64 plus a copy of one device's real and
imaginary parts during the reduction. For more detectors, compute the h-maps of groups of detectors (one
`Observation` per group) and concatenate the results.

The other buffers are bounded: pointings and pixels are transferred in blocks of samples, and in
interpolated mode the anchors of each block are computed on the host and uploaded with it (a few MB), so
their memory does not grow with the length of the observation.

## Running the GPU tests

The GPU tests are kept out of the main test suite and CI. From the repository root, with Julia 1.10 or
later:

```bash
julia --project=test/gpu -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'   # once
julia -t 4 --project=test/gpu test/gpu/runtests.jl
```

`Pkg.develop` makes the test environment use this checkout of Falcons; the test script stops if another
copy is loaded.
