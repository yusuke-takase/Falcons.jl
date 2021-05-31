## Mapmake

Mapmaking is available for computing hitmap and crosslink maps.
```julia
Mapmaking(ScanningStrategy(), split_num::Int)
```
This function splits the observation time specified by `ScanningStrategy()` by the number specified by `split_num`, calculates the pointing TOD, and then creates a map using it.
The reason for the split calculation is to avoid overloading the memory by calculating a huge amount of TOD data at once.

For example,
```julia
outmap = Mapmaking(ScanningStrategy(), 12)
```
If this is set up, TOD will be calculated every month and stored in the map.

The `outmap` is a matrix, the contents of which are as follows.
```
outmap[1]: Hitmap
outmap[2]: Crosslink map (n=1)
outmap[3]: Crosslink map (n=2)
outmap[4]: Crosslink map (n=3)
outmap[5]: Crosslink map (n=4)
```
`n` is the spin, a measure to define the scan-derived systematic effect. It is described in detail in this [paper](https://arxiv.org/abs/2008.00011).