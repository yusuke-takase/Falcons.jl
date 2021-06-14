## Scanning strategy

To define the scanning strategy for a satellite, set the `ScanningStrategy` structure.
```julia
mutable struct ScanningStrategy{T<:AbstractFloat, I<:Int, AA<:AbstractArray{T}, AS<:AbstractString}
    nside::I
    times::I
    sampling_rate::I
    alpha::T
    beta::T
    prec_period::T
    spin_rpm::T
    hwp_rpm::T
    FP_theta::AA
    FP_phi::AA
    start_point::AS
end
```

To set it up, define constructor `ss` as follows.
```julia
nside = 128
times = year #[sec]
sampling_rate = 1 #[Hz]
FP_theta = [0.0] #The angle with respect to the boresight, 0 degree represents the boresight.
FP_phi = [0.0]
alpha = 55.0 #[degree]
beta = 60.0 #[degree]
prec_period = 180.22 #[min]
spin_rpm = 0.04 #[rpm]
hwp_rpm = 0.05 #[rpm]
start_point = "pole" #You can choose "pole" or "equator"

ss = ScanningStrategy(
    nside,
    times,
    sampling_rate,
    alpha,
    beta,
    prec_period,
    spin_rpm,
    hwp_rpm,
    FP_theta,
    FP_phi,
    start_point)
```
Once the ss is set, the internal values can be accessed and modified as follows

```julia
@show ss.nside #You can see a value
ss.nside = 256 #You can change a value
```

## Generate pointing TOD
The information about the orientation of a satellite at a certain time is called pointing. The pointing is defined by $(\theta, \phi, \psi)$, where $\theta$ and $\phi$ are parameters of the 3D polar coordinates and $\psi$ is the angle between the scan direction and the meridian of the sky.

Once the scanning strategy is determined, computing the pointing is straightforward.
```julia
theta_tod, phi_tod, psi_tod, time_array = get_pointings(ss::ScanningStrategy, start::Int, stop::Int)
pix_tod, psi_tod, time_array = get_pointing_pixels(ss::ScanningStrategy, start::Int, stop::Int)
```
Enter an integer value for the time to be calculated in the `start` and `stop` fields.

`theta_tod` and `phi_tod` contain the pointing data in chronological order, and `psi_tod` contains the scan angle according to the [COSMO(HEALPix)](https://lambda.gsfc.nasa.gov/product/about/pol_convention.cfm) definition. And `time_array` contains the time used in the calculation.

The return values of these two functions are tuples, and can be combined into a single variable as an array, as shown below. In this case, the values will be stored in the following order.
```julia
pointing_TOD = get_pointings(ss::ScanningStrategy, start::Int, stop::Int)

pointing_TOD[1]: theta_tod
pointing_TOD[2]: phi_tod
pointing_TOD[3]: psi_tod
pointing_TOD[4]: time_array
```

`get_pointing_pixels()` does not allocate $\theta$ and $\phi$ arrays internally, but only allocates the minimum number of arrays needed to allocate the pixel TOD.
Therefore, it runs faster than `get_pointings()`. In fact, `get_pointing_pixels()` is executed inside `ScanningStrategy2map()`.

Now, the pointing TOD is calculated every month and stored in the map each time.
```julia
ss.FP_theta = [0.0, 10.0]
ss.FP_phi = [0.0, 90.0]
```
In this operation, the first component of the array represents the center of the focal plane, i.e. the boresight. On the other hand, the second component represents the detector that observes a point in the sky 10 degrees in theta direction and 90 degrees in phi direction away from the boresight.

If you want to build a focal plane, just substitute its configuration into these arrays and you will be able to compute multi-channel pointing.



