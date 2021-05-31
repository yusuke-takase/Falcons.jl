```@meta
CurrentModule = Falcons
```

# Falcons
This is the documentation of [Falcons](https://github.com/yusuke-takase/Falcons.jl)(Fast Algorithm for Locus Computing ON the Sky), a package for fast simulation of satellite observations.

CMB polarimetric satellites such as [LiteBIRD](http://litebird.jp/) set up an appropriate scan strategy to reduce systematic errors, and Falcons can calculate all-sky hit maps, crosslink maps, etc. by simply setting the parameters of this scan strategy. It also supports observations with multi-channel detectors.

## Installation
From the Julia REPL, run
```julia
import Pkg
Pkg.add("Falcons")
```
## Tutorial
The tutorial is available on the github page in jupyternotebook format.
Please refer to [here](https://github.com/yusuke-takase/Falcons.jl/tree/master/tutorial).

## Documentation
The documentation was built using [Documenter.jl](https://github.com/JuliaDocs).

```@autodocs
Modules = [Falcons]
```
