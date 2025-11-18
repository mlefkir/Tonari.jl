module Tonari

using Random
using FFTW
using Distributions
using Interpolations
using Unitful
using StatsBase
using ProgressMeter

include("utils.jl")
include("psd.jl")
include("simulate.jl")
include("periodogram.jl")
include("csd.jl")
include("timeseries.jl")
include("iccf.jl")

export PowerSpectralDensity,
    ContinuumPowerSpectrum,
    Model,
    PowerLaw,
    SingleBendingPowerLaw,
    DoubleBendingPowerLaw,
    PowerLaw,
    Lorentzian,
    QPO,
    evaluate,
    separate_psd, periodogram, Simulation, timmer_koenig, cross_correlate,
    time_series_sanity_checks,
    CrossSpectralDensity,
    ConstantTimeLag,
    cross_periodogram,
    ConstantPhaseLag,
    fill_gaps,
    TimeSeriesData,
    IrregularTimeStamps, RegularTimeStamps, TimeSeries, sample
end
