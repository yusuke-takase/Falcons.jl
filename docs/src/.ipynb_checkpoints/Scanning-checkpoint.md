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
ss.times = year #[sec]
ss.sampling_rate = 1 #[Hz]
ss.FP_theta = [0.0] #The angle with respect to the boresight, 0 degree represents the boresight.
ss.FP_phi = [0.0]
ss.alpha = 55 #[degree]
ss.beta = 60 #[degree]
ss.prec_period = 180.22 #[sec]
ss.spin_rpm = 0.04 #[rpm]
ss.hwp_rpm = 0.05 #[rpm]
ss.start_point = "pole" #You can choose "pole" or "equator"
```
No default values are set, so be sure to specify values for all variables by the user.

## Generate pointing TOD
The information about the orientation of a satellite at a certain time is called pointing. The pointing is defined by $(\theta, \phi, \psi)$, where $\theta$ and $\phi$ are parameters of the 3D polar coordinates and $\psi$ is the angle between the scan direction and the meridian of the sky.

Once the scanning strategy is determined, computing the pointing is straightforward.
```julia
theta_tod, phi_tod, psi_tod, time_array = get_pointings(ScanningStrategy(), start::Int, stop::Int)
pix_tod, psi_tod, time_array = get_pointing_pixels(ScanningStrategy(), start::Int, stop::Int)
```
Enter an integer value for the time to be calculated in the `start` and `stop` fields.

`theta_tod` and `phi_tod` contain the pointing data in chronological order, and `psi_tod` contains the scan angle according to the [COSMO(HEALPix)](https://lambda.gsfc.nasa.gov/product/about/pol_convention.cfm) definition. And `time_array` contains the time used in the calculation.

`get_pointing_pixels()` does not allocate $\theta$ and $\phi$ arrays internally, but only allocates the minimum number of arrays needed to allocate the pixel TOD.
Therefore, it runs faster than `get_pointings()`. In fact, `get_pointing_pixels()` is executed inside `ScanningStrategy2map()`.

Now, the pointing TOD is calculated every month and stored in the map each time.
```julia
ss.FP_theta = [0.0, 10.0]
ss.FP_phi = [0.0, 90.0]
```
In this operation, the first component of the array represents the center of the focal plane, i.e. the boresight. On the other hand, the second component represents the detector that observes a point in the sky 10 degrees in theta direction and 90 degrees in phi direction away from the boresight.

If you want to build a focal plane, just substitute its configuration into these arrays and you will be able to compute multi-channel pointing.



