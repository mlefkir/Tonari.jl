# Tonari

[![Documentation](https://github.com/mlefkir/Tonari.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/mlefkir/Tonari.jl/actions/workflows/documentation.yml) [![Build](https://github.com/mlefkir/Tonari.jl/actions/workflows/testbuild.yml/badge.svg)](https://github.com/mlefkir/Tonari.jl/actions/workflows/testbuild.yml) [![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

[![codecov](https://codecov.io/gh/mlefkir/Tonari.jl/branch/periodograms/graph/badge.svg?token=JFLB0ZE7ZY)](https://codecov.io/gh/mlefkir/Tonari.jl)

Tonari is a Julia package for the analysis of time series with an emphasis on astronomical time series. I am implementing the methods and functions I need along the way, so it is a work in progress.

Here are the current features:

- Simulation of time series with power spectral density (PSD) models (e.g. bending power law, Lorentzian, etc.)
- Simulation of bivariate time series with a cross-spectral density (CSD) model
- Periodogram computation
- Cross-periodogram for two time series with coherence, phase, and lag
- Interpolation of time series with gaps + randomisation
