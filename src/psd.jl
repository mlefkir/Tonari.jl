abstract type Model end
abstract type PowerSpectralDensity <: Model end
abstract type BendingPowerLaw <: PowerSpectralDensity end

@doc raw""" 
     SingleBendingPowerLaw(A, α₁, f₁, α₂)

Single bending power law model for the power spectral density

- `A`: the amplitude
- `α₁`: the first power law index
- `f₁`: the first bend frequency
- `α₂`: the second power law index

```math
\mathcal{P}(f) =  A \frac{(f/f₁)^{-α₁}}{1 + (f / f₁)^{α₂ - α₁}}
```

"""
struct SingleBendingPowerLaw{T<:Real} <: BendingPowerLaw
    A::T
    α₁::T
    f₁::T
    α₂::T
end

SingleBendingPowerLaw(α₁::T, f₁::T, α₂::T) where {T<:Real} = SingleBendingPowerLaw{T}(1.0, α₁, f₁, α₂)

@doc raw""" 
     DoubleBendingPowerLaw(α₁, f₁, α₂, f₂, α₃)

Double bending power law model for the power spectral density

- `α₁`: the first power law index
- `f₁`: the first bend frequency
- `α₂`: the second power law index
- `f₂`: the second bend frequency
- `α₃`: the third power law index

```math
\mathcal{P}(f) =  A\frac{(f/f₁)^{-α₁}}{1 + (f / f₁)^{α₂ - α₁}}\frac{1}{1 + (f / f₂)^{α₃ - α₂}}
```
"""
struct DoubleBendingPowerLaw{T<:Real} <: BendingPowerLaw
    A::T
    α₁::T
    f₁::T
    α₂::T
    f₂::T
    α₃::T
end

DoubleBendingPowerLaw(α₁::T, f₁::T, α₂::T, f₂::T, α₃::T) where {T<:Real} = DoubleBendingPowerLaw{T}(1.0, α₁, f₁, α₂, f₂, α₃)

@doc raw""" 
    Lorentzian(A, γ, f₀)

Lorentzian model for the power spectral density

- `A`: the amplitude
- `γ`: the width of the peak
- `f₀`: the central frequency

```math
\mathcal{P}(f) =  \frac{A}{4\pi^2 (f - f₀)^2 + γ^2}
```
"""
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
    return psd.A * (f / psd.f₁)^(-psd.α₁) / (1 + (f / psd.f₁)^(psd.α₂ - psd.α₁)) / (1 + (f / psd.f₂)^(psd.α₃ - psd.α₂))
end

function calculate(f, psd::SingleBendingPowerLaw)
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
