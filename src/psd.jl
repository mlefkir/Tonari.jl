abstract type Model end
abstract type PowerSpectralDensity <: Model end
abstract type BendingPowerLaw <: PowerSpectralDensity end


struct SimpleBendingPowerLaw{T<:Real} <: BendingPowerLaw
    A::T
    α₁::T
    f₁::T
    α₂::T
end

SimpleBendingPowerLaw(α₁::T, f₁::T, α₂::T) where {T<:Real} = SimpleBendingPowerLaw{T}(1.0, α₁, f₁, α₂)

struct DoubleBendingPowerLaw{T<:Real} <: BendingPowerLaw
    A::T
    α₁::T
    f₁::T
    α₂::T
    Δf::T
    α₃::T
end

SimpleBendingPowerLaw(α₁::T, f₁::T, α₂::T, Δf::T, α₃::T) where {T<:Real} = SimpleBendingPowerLaw{T}(1.0, α₁, f₁, α₂, Δf, α₃)

struct Lorentzian{TA<:Real,Tγ<:Real,Tf₀<:Real} <: PowerSpectralDensity
    A::TA
    γ::Tγ
    f₀::Tf₀
end

Lorentzian(γ::Tγ, f₀::Tf₀) where {Tγ<:Real,Tf₀<:Real} = Lorentzian(1.0, γ, f₀)


function calculate(f, psd::Lorentzian)
    return psd.A ./ (4π^2 .* (f .- psd.f₀) .^ 2 + psd.γ^2)
end

function calculate(f, psd::DoubleBendingPowerLaw)
    return psd.A * (f / psd.f₁)^(-psd.α₁) / (1 + (f / psd.f₁)^(psd.α₂ - psd.α₁)) / (1 + (f / (psd.f₁ * psd.Δf))^(psd.α₃ - psd.α₂))
end

function calculate(f, psd::SimpleBendingPowerLaw)
    return psd.A * (f / psd.f₁)^(-psd.α₁) / (1 + (f / psd.f₁)^(psd.α₂ - psd.α₁))
end


struct SumOfPowerSpectralDensity{Tp<:Vector{<:PowerSpectralDensity}} <: PowerSpectralDensity
    psd::Tp
end

function Base.:+(a::PowerSpectralDensity, b::PowerSpectralDensity)
    SumOfPowerSpectralDensity([a, b])
end

function calculate(f, model::SumOfPowerSpectralDensity)
    sum(p(f) for p in model.psd)
end

(psd::PowerSpectralDensity)(f) = calculate.(f, Ref(psd))
