# Falcons
<img src="https://user-images.githubusercontent.com/83496454/132131349-2c56c95c-ab8c-45ec-889a-1fe8e68d4d83.png" width="600">

This logo was created with the help of our dear collaborator Jonathan Aumont and Léo Vacher. I would like to thank you from the bottom of my heart.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yusuke-takase.github.io/Falcons.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yusuke-takase.github.io/Falcons.jl/dev)
[![Build Status](https://travis-ci.com/yusuke-takase/Falcons.jl.svg?branch=master)](https://travis-ci.com/yusuke-takase/Falcons.jl)

Welcome to Falcons(Fast Algorithm for Locus Computing ON the Sky).
Falcons is a software that rapidly calculates the TOD(Time Ordered Data) of pointing information required for satellite observation simulations.

It supports multi-detectors and can construct a focal plane by specifying (theta,phi) arrays centered on the boresight.
The mapmaking function creates a hitmap and crosslink map from the rapidly obtained time series pointing information.

Falcons is fast enough to run on a laptop, but can be easily used on a supercomputer. 
Since the user can specify the amount of memory to occupy according to the available memory on the node, it is very suitable for jobs that are submitted in large quantities.

![Figure](https://user-images.githubusercontent.com/83496454/119337906-532ff680-bcca-11eb-9b8c-bde7a376c6e6.gif)

This is how the sky is scanned by 607 detectors computed by Falcons.

## Installation
From the Julia REPL, run

```julia
import Pkg
Pkg.add("Falcons")
```

## Threading 
Falcons uses multithreading technology internally when computing pointing TOD. It uses the `@threads` macro provided by Julia, and to use this feature, users need to make the following declarations in the configuration file of the computing environment (ex. `.bashrc`) in advance.
```bash
export JULIA_NUM_THREADS=4
```
In this example, you will use 4 threads. You can change the number of threads according to your computing environment.

## Usage example
See the tutorial [here](https://github.com/yusuke-takase/Falcons.jl/tree/master/tutorial) for details on how to use it.
Refer to the [documentation](https://yusuke-takase.github.io/Falcons.jl/dev/) for more examples.

