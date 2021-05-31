## Scanning strategy

To define the scanning strategy for a satellite, set the `ScanningStrategy` structure.
```julia
mutable struct ScanningStrategy
    nside::Int
    times::Int
    sampling_rate::Int
    alpha::AbstractFloat
    beta::AbstractFloat
    prec_period::AbstractFloat
    spin_rpm::AbstractFloat
    hwp_rpm::AbstractFloat
    FP_theta::AbstractArray{AbstractFloat,1}
    FP_phi::AbstractArray{AbstractFloat,1}
    start_point::AbstractString
    ScanningStrategy() = new()
end
```

To set it up, define constructor `ss` as follows.
```julia
ss = ScanningStrategy()
```
You can assign a value to `ss` by accessing it as follows.
```julia
ss.nside = 128
ss.times = 100
```
No default values are set, so be sure to specify values for all variables by the user.

## Generate pointing TOD

Once the scanning strategy is determined, computing the pointing is straightforward.
```julia
theta_tod, phi_tod, psi_tod = get_scan_tod(ScanningStrategy(), start::Int, stop::Int)
```
Enter an integer value for the time to be calculated in the `start` and `stop` fields.

`theta_tod` and `phi_tod` contain the pointing data in chronological order, and `psi_tod` contains the scan angle according to the IAU definition.
