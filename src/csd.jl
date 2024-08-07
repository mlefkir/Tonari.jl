abstract type PhaseModel <: Model end

@doc raw""" 
    CrossSpectralDensity(𝓟₁::PowerSpectralDensity, 𝓟₂::PowerSpectralDensity, Δφ::PhaseModel)

Cross spectral density model, stores the power spectral density models of the individual processes
and the phase lag between them
"""
struct  CrossSpectralDensity <: Model
    𝓟₁::PowerSpectralDensity
    𝓟₂::PowerSpectralDensity
    Δφ::PhaseModel
end

@doc raw"""
    ConstantTimeLag

Constant time lag model, stores the time lag between the two processes
"""
struct ConstantTimeLag <: PhaseModel
    Δt::Real
end

function calculate(Δφ::ConstantTimeLag, f)
    return Δφ.Δt
end
