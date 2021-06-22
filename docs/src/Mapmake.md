## Mapmake

Mapmaking is available for computing hitmap and crosslink maps.
```julia
ScanningStrategy2map(ss::ScanningStrategy, division::Int)
```
This function splits the observation time specified by `ScanningStrategy` by the number specified by `division`, calculates the pointing TOD, and then creates a map using it.
The reason for the split calculation is to avoid overloading the memory by calculating a huge amount of pointing TOD data at once.
Inside this function, `get_pointing_pixels()` is being executed.

For example,
```julia
outmap = ScanningStrategy2map(ss::ScanningStrategy, 12)
```
Now, the pointing TOD is calculated every month and stored in the map each time.

The `outmap` is a matrix, the contents of which are as follows.
```
outmap[1]: Hitmap
outmap[2]: Crosslink map (n=1)
outmap[3]: Crosslink map (n=2)
outmap[4]: Crosslink map (n=3)
outmap[5]: Crosslink map (n=4)
```
The crosslink is represented by the following equation

$Crosslink(n)=<\sin(n\psi_i)>^2+<\cos(n\psi_i)>^2$.
The subscript $i$ indicates that this is the $i$-th time series data within a certain pixel, and $< >$ indicates that the average should be taken.

The $n$ is the spin, a measure to define the scan-derived systematic effect. It is described in detail in this [paper](https://arxiv.org/abs/2008.00011).

## Crosslink with half-wave plate
The next generation of CMB polarimetry satellites, such as LiteBIRD, will be equipped with a half-wave plate(HWP). Falcons is also possible to create a crosslink map that takes into account the optical axis of the HWP.
The crosslink with HWP is calculated by the following equation
$Crosslink_{hwp}(n)=<\sin(4\rho-n\psi)>^2+<\cos(4\rho-n\psi)>^2$.