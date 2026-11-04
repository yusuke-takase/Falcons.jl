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

CUDA is a weak dependency: Falcons works without it, and calling a GPU backend without `using CUDA` raises
an `ArgumentError`.

## Running the GPU tests

The GPU tests are kept out of the main test suite and CI:

```bash
julia --project=test/gpu -e 'using Pkg; Pkg.instantiate()'
julia -t 4 --project=test/gpu test/gpu/runtests.jl
```
