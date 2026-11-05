# Falcons

<img src="https://user-images.githubusercontent.com/83496454/132532967-9c2f0e19-d920-4b94-863f-93236e093ff9.png" width="600">

This logo was created with the help of our dear collaborator Jonathan Aumont.

![Build Status](https://github.com/yusuke-takase/Falcons.jl/actions/workflows/test.yml/badge.svg?branch=master)

Welcome to Falcons (Framework for Astrophysical Locus Computing ON the Sky).
Falcons rapidly computes the pointings of spinning CMB satellites and the maps derived from them: hit maps
and the spin moments of the scan (h-maps, also called cross-link maps), for any number of detectors placed
on the focal plane.

Falcons follows the conventions and names of [litebird_sim](https://github.com/litebird/litebird_sim): its
tests check that both codes give the same pointings, HWP angles, hit maps and h-maps, it reads
litebird_sim's instrument databases (IMO) and writes h-maps in litebird_sim's file layout. It runs
multithreaded on a CPU and, optionally, on one or several CUDA GPUs.

![Figure](https://user-images.githubusercontent.com/83496454/155742440-294f6b97-1305-43ac-8d57-8534eeab7005.gif)

## Installation

```julia
import Pkg
Pkg.add(url = "https://github.com/yusuke-takase/Falcons.jl")
```

Falcons uses `Threads.@threads`: start Julia with several threads, e.g. `julia -t 8` (or set
`JULIA_NUM_THREADS`). For the GPU backend, also install [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) and
run `using CUDA` before calling Falcons with `backend = GPUBackend()`.

## Usage example

```julia
using Falcons

scan, inst = satellite_preset(:LiteBIRD)
dets = [
    DetectorInfo(name = "boresight", sampling_rate_hz = 19.0),
    DetectorInfo(name = "offset", sampling_rate_hz = 19.0, quat = detector_quat(deg2rad(5), deg2rad(30))),
]
obs = Observation(scan, inst, dets; start_time_s = 0.0, duration_s = 86400.0,
                  hwp = IdealHWP(ang_speed_radpsec = 46 * 2π / 60))

pointings, hwp_angle = get_pointings(obs)        # (θ, φ, ψ) of every sample and detector
hits = get_hitmap(obs; nside = 128)
res = make_h_maps(obs; nside = 128, n_m_couples = nm_grid(0:4, (-4, 0, 4)))
save_h_maps(res, "h_maps")
```

See the [tutorial](tutorial/falcons_tutorial.ipynb) and the documentation for the conventions, the IMO
interface, time-dependent attitudes and the GPU backend. To build the documentation locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path = pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl     # output in docs/build
```

## How to cite the code

If you use this software, please cite [Y. Takase et al, 2024, JCAP](https://iopscience.iop.org/article/10.1088/1475-7516/2024/12/036) for your publication.
