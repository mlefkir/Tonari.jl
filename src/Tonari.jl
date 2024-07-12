module Tonari
using Random
using FFTW
using Distributions

include("psd.jl")
include("simulate.jl")
include("periodogram.jl")

export periodogram, Simulation, PowerSpectralDensity, SingleBendingPowerLaw, DoubleBendingPowerLaw, Lorentzian, timmer_koenig, calculate, sample

end
