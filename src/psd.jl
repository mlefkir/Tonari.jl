abstract type Model end
abstract type PowerSpectralDensity <: Model end
abstract type BendingPowerLaw <: PowerSpectralDensity end


struct SimpleBendingPowerLaw{T<:Real} <: BendingPowerLaw
    α₁::T
    f₁::T
    α₂::T
end

struct DoubleBendingPowerLaw{T<:Real} <: BendingPowerLaw
    α₁::T
    f₁::T
    α₂::T
    Δf::T
    α₃::T
end


struct Lorentzian{T<:Real} <: PowerSpectralDensity
    A::T
    γ::T
    f₀::T
end


function calculate(f, psd::Lorentzian)
    return psd.A ./ (1 .+ (f .- psd.f₀) .^ 2 / psd.γ^2)

end


function calculate(f, psd::DoubleBendingPowerLaw)
    return (f / psd.f₁)^(-psd.α₁) / (1 + (f / psd.f₁)^(psd.α₂ - psd.α₁)) / (1 + (f / (psd.f₁ * psd.Δf))^(psd.α₃ - psd.α₂))

end

function calculate(f, psd::SimpleBendingPowerLaw)
    return (f / psd.f₁)^(-psd.α₁) / (1 + (f / psd.f₁)^(psd.α₂ - psd.α₁))
end
