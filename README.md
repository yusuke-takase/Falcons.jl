# Falcons
Welcome to FLCONS(Fast Algorithm for Locus Computing ON the Sky).
FALCONS is a software that rapidly calculates the tod of pointing information required for satellite observation simulations. FALCONS supports multi-detectors and can construct a focal plane by specifying (theta,phi) arrays centered on the boresight.
The mapmaking function creates a crosslink map from the rapidly obtained time series pointing information.

![Multi-detector's trajectory](https://user-images.githubusercontent.com/83496454/119337906-532ff680-bcca-11eb-9b8c-bde7a376c6e6.gif)

## Installation
From the Julia REPL, run

```julia
import Pkg
Pkg.add("Falcons.jl")
```

## Usage example
The program in jupyter format is uploaded to tutorial.
Please refer to it.


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yusuke-takase.github.io/Falcons.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yusuke-takase.github.io/Falcons.jl/dev)
[![Build Status](https://travis-ci.com/yusuke-takase/Falcons.jl.svg?branch=master)](https://travis-ci.com/yusuke-takase/Falcons.jl)
