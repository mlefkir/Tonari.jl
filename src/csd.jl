abstract type PhaseModel <: Model end

@doc raw"""
	CrossSpectralDensity(𝓟₁::PowerSpectralDensity, 𝓟₂::PowerSpectralDensity, Δφ::PhaseModel)

Cross spectral density model, stores the power spectral density models of the individual processes
and the phase lag between them.
"""
struct CrossSpectralDensity <: Model
    𝓟₁::PowerSpectralDensity
    𝓟₂::PowerSpectralDensity
    Δφ::PhaseModel
end

@doc raw"""
	ConstantTimeLag

Constant time lag model, stores the time lag between the two processes
"""
struct ConstantTimeLag <: PhaseModel
    Δτ::Real
end


function evaluate(Δφ::ConstantTimeLag, f)
    return Δφ.Δτ
end

@doc raw"""
	ConstantPhaseLag

Constant phase lag model, stores the phase lag between the two processes

# Fields
- `τ₀::Real`: Time lag at f₀
- `f₀::Real`: Reference frequency
"""
struct ConstantPhaseLag <: PhaseModel
    τ₀::Real
    f₀::Real
end

function evaluate(Δφ::ConstantPhaseLag, f)
    return @. Δφ.τ₀ * Δφ.f₀ / f
end


(Δφ::PhaseModel)(f) = calculate.(f, Ref(Δφ))
