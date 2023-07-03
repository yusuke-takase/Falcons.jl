# Falcons
<img src="https://user-images.githubusercontent.com/83496454/132532967-9c2f0e19-d920-4b94-863f-93236e093ff9.png" width="600">

This logo was created with the help of our dear collaborator Jonathan Aumont. I would like to thank you from the bottom of my heart.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yusuke-takase.github.io/Falcons.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yusuke-takase.github.io/Falcons.jl/dev)
[![Build Status](https://travis-ci.com/yusuke-takase/Falcons.jl.svg?branch=master)](https://travis-ci.com/yusuke-takase/Falcons.jl)

Welcome to Falcons(Framework for Astrophysical Locus Computing ON the Sky).
Falcons is a software that rapidly calculates the TOD(Time Ordered Data) of pointing information required for satellite observation simulations.

It supports multi-detectors and can construct a focal plane by specifying (theta,phi) arrays centered on the boresight.
The mapmaking function creates a hitmap and crosslink map from the rapidly obtained time series pointing information.

Falcons is fast enough to run on a laptop, but can be easily used on a supercomputer.
Since the user can specify the amount of memory to occupy according to the available memory on the node, it is very suitable for jobs that are submitted in large quantities.

![Figure](https://user-images.githubusercontent.com/83496454/155742440-294f6b97-1305-43ac-8d57-8534eeab7005.gif)

## Installation
From the Julia REPL, in order to install Falcons from Julia general repository, you can run

```julia
import Pkg
Pkg.add("Falcons")
```
Or you can install the Falcons directory from this GitHub by running
```julia
Pkg> add https://github.com/yusuke-takase/Falcons.jl.git
```
Then, master version is going to be install in your environment.

## Threading 
Falcons uses multithreading technology internally when computing pointing TOD. It uses the `@threads` macro provided by Julia, and to use this feature, users need to make the following declarations in the configuration file of the computing environment (ex. `.bashrc`) in advance.
```bash
export JULIA_NUM_THREADS=4
```
In this example, you will use 4 threads. You can change the number of threads according to your computing environment.

## Usage example
See the tutorial [here](https://github.com/yusuke-takase/Falcons.jl/tree/master/tutorial) for details on how to use it.
Refer to the [documentation](https://yusuke-takase.github.io/Falcons.jl/dev/) for more examples.

