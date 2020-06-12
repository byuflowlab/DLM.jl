# DLM.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/DLM.jl)
![](https://github.com/byuflowlab/DLM.jl/workflows/Run%20tests/badge.svg)

*Pure Julia implementation of the Doublet Lattice Method*

Author: Taylor McDonnell

**DLM.jl** is a pure Julia implementation of the Doublet Lattice Method.

## Package Features
- Capable of modeling multiple arbitrarily defined trapezoidal geometries
- Highly accurate 12-term exponential approximation of the kernel integrals
- Multiple options for fits of the linear and nonlinear terms of the kernel function
 - With the parabolic fit:
  - Capable of accurately handling panel aspect ratios of up to 3
 - With the quartic fit:
  - Capable of accurately handling panel aspect ratios of up to 6-10
- Simple result visualization with [WriteVTK](https://github.com/jipolanco/WriteVTK.jl)
- Validated against published Doublet Lattice Method results.
- Allocation-free code

## Installation

Enter the package manager by typing `]` and then run the following:

```julia
pkg> add https://github.com/byuflowlab/DLM.jl
```

## Performance

This package was originally created based on the Python/Fortran DLM implementation in [dlm4py](https://github.com/gjkennedy/dlm4py). The DLM implemented in this package is roughly 50% slower than that Fortran implementation.  However, since this package is implemented in pure Julia, it can be used with automatic differentiation as implemented by the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) or [ReverseDiff](https://github.com/JuliaDiff/ReverseDiff.jl) packages, or used with other dynamic language features.

## Usage

See the [example](@ref Example)
