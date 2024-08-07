module Tonari
using Random
using FFTW
using Distributions

include("psd.jl")
include("simulate.jl")
include("periodogram.jl")
include("csd.jl")

export periodogram, Simulation, PowerSpectralDensity, SingleBendingPowerLaw, DoubleBendingPowerLaw, Lorentzian, timmer_koenig, calculate, sample, CrossSpectralDensity, ConstantTimeLag, cross_periodogram
end
