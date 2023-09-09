## Scanning strategy

To define the scanning strategy for a satellite, set the `ScanningStrategy` structure.
```julia
mutable struct ScanningStrategy{T<:AbstractFloat, I<:Int, AS<:AbstractString}
    nside::I
    duration::I
    sampling_rate::T
    alpha::T
    beta::T
    gamma::T
    prec_rpm::T
    spin_rpm::T
    hwp_rpm::T
    start_point::AS
    start_angle::T
    coord::AS
    quat::Vector{Vector{Float64}}
    name::Vector{String}
    info::DataFrame
end
```
You can generate `ScanningStrategy` structure using `gen_ScanningStrategy()` function.


This initial value of component of `ScanningStrategy` can be changed by specifying the `gen_ScanningStrategy()` argument when declaring `ss` like below.
```julia
ss = gen_ScanningStrategy(alpha=60, prec_rpm=0.001, sampling_rate=5.0)
```
In addition, you can directly access ss and interactively change its value.
```julia
ss.nside = 256 #You can change a value
ss.spin_rpm = 0.04
```

## Generate pointings
The information about the orientation of a satellite at a certain time is called pointing. The pointing is defined by $(\theta, \phi, \psi)$, where $\theta$ and $\phi$ are parameters of the 3D polar coordinates and $\psi$ is the angle between the scan direction and the meridian of the sky.

Once the scanning strategy is determined, computing the pointing is straightforward.
```julia
pointings = get_pointings(ss::ScanningStrategy, start::Int, stop::Int)
pix_tod, psi_tod, time_array = get_pointing_pixels(ss::ScanningStrategy, start::Int, stop::Int)
```
Enter an integer value for the time to be calculated in the `start` and `stop` fields.
The `get_pointings()` returns $(\theta, \phi, \psi)$ that indicate their pointing in chronological order. These are stored as an array of dictionary type, and can be accessed as follows.
```julia
In[]: pointings["theta"]
Out[]: 14400×1 Matrix{Float64}:
 1.1344640137963142
 1.134461356131033
 1.1344533831585348
 1.1344400949488544
...
```

The `pointings["psi"]` contains the scanning angle according to the [COSMO(HEALPix)](https://lambda.gsfc.nasa.gov/product/about/pol_convention.cfm) definition. And `pointings["time"]` contains the time used in the calculation.

The `get_pointing_pixels()` does not allocate $\theta$ and $\phi$ arrays internally, but only allocates the minimum number of arrays needed to allocate the pixel TOD.
Therefore, it runs faster than `get_pointings()`. In fact, `get_pointing_pixels()` is executed inside `ScanningStrategy2map()`.

Now, the pointing TOD is calculated every month and stored in the map each time.
```julia
ss.FP_theta = [0.0, 10.0]
ss.FP_phi = [0.0, 90.0]
```
In this operation, the first component of the array represents the center of the focal plane, i.e. the boresight. On the other hand, the second component represents the detector that observes a point in the sky 10 degrees in theta direction and 90 degrees in phi direction away from the boresight.

If you want to build a focal plane, just substitute its configuration into these arrays and you will be able to compute multi-channel pointings.


