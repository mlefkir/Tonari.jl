# Tonari

[![Documentation](https://github.com/mlefkir/Tonari.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/mlefkir/Tonari.jl/actions/workflows/documentation.yml) [![Build](https://github.com/mlefkir/Tonari.jl/actions/workflows/testbuild.yml/badge.svg)](https://github.com/mlefkir/Tonari.jl/actions/workflows/testbuild.yml) [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

[![codecov](https://codecov.io/gh/mlefkir/Tonari.jl/branch/periodograms/graph/badge.svg?token=JFLB0ZE7ZY)](https://codecov.io/gh/mlefkir/Tonari.jl)

Tonari is a Julia package for the analysis of time series with an emphasis on astronomical time series. I am implementing the methods and functions I need along the way, so it is a work in progress.

## Installation

Since it's Julia package, you need... Julia! See the official documentation here: [https://julialang.org/](https://julialang.org/).

To enter Pkg mode by type `]` in the julia REPL
```julia
pkg>
```

To install Tonari type:
```julia
pkg> add Tonari
```

and then load it with:
```julia
using Tonari
```

## Documentation

The documentation is available here: [https://mlefkir.github.io/Tonari.jl/](https://mlefkir.github.io/Tonari.jl/)

## Current features

- Power spectral density models: power-law with zero, one or two bends, Lorentzian...
- Simulation of time series with power spectral density (PSD) models (e.g. bending power law, Lorentzian, etc.)
- Simulation of bivariate time series with a cross-spectral density (CSD) model
- Periodogram computation
- Cross-periodogram for two time series with coherence, phase, and lag
- Interpolation of time series with gaps + randomisation
- Delay estimation with the interpolated cross-correlation function (ICCF)
- Time series structure to represent regularly and irregularly sampled time series

## Ideas/possible features

### for the 0.3.x
- Periodogram fitting with Whittle likelihood, essentially connect sampler and a good likelihood definition.
- Extracting light curves from event lists, mainly from X-ray observations.
### for the 0.4.x and later
- Operation on light curves, concatenate light curves with dates, i.e. making a big light curve.
- Bispectrum computation for fun.