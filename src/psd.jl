abstract type Model end
abstract type PowerSpectralDensity <: Model end
abstract type BendingPowerLaw <: PowerSpectralDensity end

@doc raw"""
	QPO(S₀, f₀,A Q)

QPO model

- `S₀`: the amplitude at the peak
- `f₀`: the central frequency
- `Q`: quality factor

```math
\mathcal{P}(f) =  \frac{S_0 {f_0}^4 } {\left(f^ 2 -{f_0}^2\right)^ 2 + f^2 {f_0}^2 /  Q^2 }
```
"""
struct QPO{TS₀ <: Real, Tf₀ <: Real, TQ <: Real} <: PowerSpectralDensity
    S₀::TS₀
    f₀::Tf₀
    Q::TQ

end

QPO(f₀::Tf₀, Q::TQ) where {Tf₀ <: Real, TQ <: Real} = QPO(1.0, f₀, Q)

@doc raw"""
     PowerLaw(α)

Power law model for the power spectral density

- `α`: the power law index

```math
\mathcal{P}(f) = A f^{-α}
```

"""
struct PowerLaw{T <: Real} <: BendingPowerLaw
    A::T
    α::T
end

PowerLaw(α::T) where {T <: Real} = PowerLaw{T}(1.0, α)


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
struct SingleBendingPowerLaw{T <: Real} <: BendingPowerLaw
    A::T
    α₁::T
    f₁::T
    α₂::T
end

SingleBendingPowerLaw(α₁::T, f₁::T, α₂::T) where {T <: Real} = SingleBendingPowerLaw{T}(1.0, α₁, f₁, α₂)

@doc raw"""
     DoubleBendingPowerLaw(A, α₁, f₁, α₂, f₂, α₃)

Double bending power law model for the power spectral density

- `A` : the amplitude
- `α₁`: the first power law index
- `f₁`: the first bend frequency
- `α₂`: the second power law index
- `f₂`: the second bend frequency
- `α₃`: the third power law index

```math
\mathcal{P}(f) =  A\frac{(f/f₁)^{-α₁}}{1 + (f / f₁)^{α₂ - α₁}}\frac{1}{1 + (f / f₂)^{α₃ - α₂}}
```
"""
struct DoubleBendingPowerLaw{T <: Real} <: BendingPowerLaw
    A::T
    α₁::T
    f₁::T
    α₂::T
    f₂::T
    α₃::T
end

DoubleBendingPowerLaw(α₁::T, f₁::T, α₂::T, f₂::T, α₃::T) where {T <: Real} = DoubleBendingPowerLaw{T}(1.0, α₁, f₁, α₂, f₂, α₃)

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
struct Lorentzian{TA <: Real, Tγ <: Real, Tf₀ <: Real} <: PowerSpectralDensity
    A::TA
    γ::Tγ
    f₀::Tf₀
end

struct SumOfPowerSpectralDensity{Tp <: Vector{<:PowerSpectralDensity}} <: PowerSpectralDensity
    psd::Tp
end


Lorentzian(γ::Tγ, f₀::Tf₀) where {Tγ <: Real, Tf₀ <: Real} = Lorentzian(1.0, γ, f₀)

function evaluate(psd::QPO, f)
    return psd.S₀ * psd.f₀^4 ./ ((f .^ 2 .- psd.f₀ .^ 2) .^ 2 .+ f .^ 2 .* psd.f₀^2 / psd.Q^2)
end

function evaluate(psd::PowerLaw, f)
    return psd.A * (f)^(-psd.α)
end

function evaluate(psd::Lorentzian, f)
    return psd.A ./ (4π^2 .* (f .- psd.f₀) .^ 2 + psd.γ^2)
end

function evaluate(psd::DoubleBendingPowerLaw, f)
    return psd.A * (f / psd.f₁)^(-psd.α₁) / (1 + (f / psd.f₁)^(psd.α₂ - psd.α₁)) / (1 + (f / psd.f₂)^(psd.α₃ - psd.α₂))
end

function evaluate(psd::SingleBendingPowerLaw, f)
    return psd.A * (f / psd.f₁)^(-psd.α₁) / (1 + (f / psd.f₁)^(psd.α₂ - psd.α₁))
end

function evaluate(model::SumOfPowerSpectralDensity, f)
    return sum(p(f) for p in model.psd)
end

function Base.:+(a::PowerSpectralDensity, b::PowerSpectralDensity)
    return SumOfPowerSpectralDensity([a, b])
end

function Base.:+(a::PowerSpectralDensity, b::SumOfPowerSpectralDensity)
    return SumOfPowerSpectralDensity([a; b.psd])
end

function Base.:+(a::SumOfPowerSpectralDensity, b::PowerSpectralDensity)
    return SumOfPowerSpectralDensity([a.psd; b])
end

function Base.:+(a::SumOfPowerSpectralDensity, b::SumOfPowerSpectralDensity)
    return SumOfPowerSpectralDensity([a.psd; b.psd])
end

(psd::PowerSpectralDensity)(f) = evaluate.(Ref(psd), f)

"""
	 separate_psd(psd::PowerSpectralDensity)

Separate the PSD into its BendingPowerLaw components and other components if it is a sum of PSDs

# Arguments
- `psd::PowerSpectralDensity`: power spectral density or sum of PowerSpectralDensity objects

# Return
- `psd_continuum::Union{SumOfPowerSpectralDensity,PowerSpectralDensity,nothing}`: continuum part of the psd
- `psd_line::Union{PowerSpectralDensity,nothing,Vector{PowerSpectralDensity}}`: non-continuum part of the psd
"""
function separate_psd(psd::PowerSpectralDensity)
    if isa(psd, BendingPowerLaw)
        return psd, nothing
    elseif isa(psd, SumOfPowerSpectralDensity)
        cont = isa.(psd.psd, BendingPowerLaw)
        # if it's a sum of features only
        if all(cont .== false)
            return nothing, psd.psd
            # if it's a sum of bending power law only
        elseif all(cont .== true)
            return psd, nothing
        else
            if length(psd.psd[cont]) < 2
                psd_continuum = psd.psd[cont][1]
            else
                psd_continuum = SumOfPowerSpectralDensity(psd.psd[cont])
            end
            psd_line = psd.psd[.!cont]
        end
        return psd_continuum, psd_line

    else
        return nothing, psd
    end
end
