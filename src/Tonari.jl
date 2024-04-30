module Tonari
using Random
using FFTW
using Distributions

include("psd.jl")
include("simulate.jl")

export Simulation, PowerSpectralDensity, SingleBendingPowerLaw, DoubleBendingPowerLaw, Lorentzian, timmer_koenig, calculate, sample

end
