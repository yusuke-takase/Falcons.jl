# Falcons

[Falcons](https://github.com/yusuke-takase/Falcons.jl) (Framework for Astrophysical Locus Computing ON the Sky)
computes the pointings of spinning CMB satellites and the quantities derived from them: hit maps and the
spin moments of the scan ("h-maps", also called cross-link maps). It is fast enough for a laptop, runs
multithreaded on a CPU and, optionally, on one or several CUDA GPUs.

Falcons follows the conventions, type names and parameter names of
[litebird_sim](https://github.com/litebird/litebird_sim), and its tests check that both codes give the same
pointings, HWP angles, hit maps and h-maps.

## Installation

```julia
import Pkg
Pkg.add(url = "https://github.com/yusuke-takase/Falcons.jl")
```

Falcons uses `Threads.@threads`; start Julia with several threads, e.g. `julia -t 8`. For the GPU backend,
install [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) as well and run `using CUDA` (see [GPU backend](@ref)).

## Quick start

```julia
using Falcons

scan, inst = satellite_preset(:LiteBIRD)               # spin/precession and spin-boresight angle
dets = [
    DetectorInfo(name = "boresight", sampling_rate_hz = 19.0),
    DetectorInfo(name = "offset", sampling_rate_hz = 19.0, quat = detector_quat(deg2rad(5), deg2rad(30))),
]
obs = Observation(scan, inst, dets; start_time_s = 0.0, duration_s = 86400.0,
                  hwp = IdealHWP(ang_speed_radpsec = 46 * 2π / 60))

pointings, hwp_angle = get_pointings(obs)               # (3, n_samples, n_detectors): θ, φ, ψ
hits = get_hitmap(obs; nside = 128)                     # Galactic, RING
res = make_h_maps(obs; nside = 128, n_m_couples = nm_grid(0:4, (-4, 0, 4)))
h_map(res, "boresight", 2, 0)                           # h_{2,0} of one detector
save_h_maps(res, "h_maps")                              # litebird_sim's file layout
```

## Contents

- [Conventions](@ref): rotation chain, angles, HWP and h-map definitions.
- [Time-dependent rotations and external attitude](@ref): `RotQuaternion`, measured attitudes, perturbations.
- [Instrument database and parameter files](@ref): reading litebird_sim's IMO, TOML pipelines.
- [GPU backend](@ref): computing on CUDA GPUs.
- [Migrating from Falcons 1.0.1](@ref): old names and their replacements.
- [API reference](@ref).

## Citing

If you use Falcons, please cite [Y. Takase et al., JCAP 2024](https://iopscience.iop.org/article/10.1088/1475-7516/2024/12/036).
