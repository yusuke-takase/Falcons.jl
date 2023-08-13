## Mapmake
# Pointings
To obtain TOD from Sky map, it is recommended to create a Pointing structure. The elemets are defined as:
```julia
mutable struct Pointings
    x::AbstractFloat
    y::AbstractFloat
    z::AbstractFloat
    θ::AbstractFloat # position of sky
    φ::AbstractFloat # position of sky
    Ω::Int           # Sky pixel index
    ψ::AbstractFloat # Crossing angle
    ϕ::AbstractFloat # HWP angle
    ξ::AbstractFloat # mod2pi(2ϕ) + ψ
end
```
It is convenient to use the `pointings` for the actual generation.
```julia
pointings(resol::Resolution, θ, φ, ψ, ϕ)
```
# Available bolometric equations
```julia
true_signal(p::Pointings, maps::PolarizedHealpixMap, pixbuf, weightbuf)
true_signal(p::Pointings, maps::PolarizedHealpixMap)
interp_signal(p::Pointings, maps::PolarizedHealpixMap)
interp_signal(p::Pointings, maps::PolarizedHealpixMap, pixbuf, weightbuf)
```
`true_signal` returns TOD which is defined by
$d_j = I(\Omega) + Q(\Omega)\cos(2\xi_j) + U(\Omega)\sin(2\xi_j)$,
where $\Omega$ shows sky pixel and $\xi\equiv\2phi+\psi$ shows effective crossing angle which takes in to account the satellite crossing angle($\psi$) and its HWP angle($\phi$).

`interp_signal` returns interpolated TOD which is defined pointings i.e. $\theta$ and $\varphi$.

# Map-maker
`binned_mapmake` is available. 
```julia
binned_mapmake(ss::ScanningStrategy, division::Int, inputinfo::Falcons.InputInfo, signal)
```
`signal` argument requires a function which has `Pointings` and `PolarizedHealpixMap` as arguments in the order. For example,
```julia
binned_mapmake(ss::ScanningStrategy, division::Int, inputinfo::Falcons.InputInfo, true_signal)
```