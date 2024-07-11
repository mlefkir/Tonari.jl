"""
    Simulation

A struct that contains the information of a simulation of a stochastic process.

# Fields
- `model::PowerSpectralDensity`: The power spectral density of the process.
- `T::Real`: The duration of the simulated time series.
- `Δt::Real`: The sampling period, or minimum time difference between samples.
- `S_high::Real`: The factor by which the maximum frequency is multiplied for the simulation.
- `S_low::Real`: The factor by which the minimum frequency is divided for the simulation.
- `t::AbstractVector{Real}`: The time vector, in this case, the time at which the process is sampled. It is assumed that the time vector is sorted.

# Constructors
- `Simulation(model::PowerSpectralDensity, T::Real, Δt::Real, S_high::Real, S_low::Real, t::AbstractVector{Real})`: Constructs a simulation with regular sampling.
- `Simulation(model::PowerSpectralDensity, T::Real, Δt::Real)`: Constructs a simulation with regular sampling, and sets `S_high` and `S_low` to 10.0.
- `Simulation(model::PowerSpectralDensity, T::Real, Δt::Real, S_high::Real, S_low::Real)`: Constructs a simulation with the given parameters.
- `Simulation(model::PowerSpectralDensity, t::AbstractVector{Real}, S_high::Real, S_low::Real)`: Constructs a simulation with sampling pattern given by `t`.
- `Simulation(model::PowerSpectralDensity, t::AbstractVector{Real})`: Constructs a simulation with sampling pattern given by `t`, and sets `S_high` and `S_low` to 10.0.
"""

struct Simulation{Tp<:PowerSpectralDensity,Tt<:Real,Ts<:Real,Tsh<:Real,Tsl<:Real}
    model::Tp
    T::Tt
    Δt::Ts
    S_high::Tsh
    S_low::Tsl
    t::AbstractVector{Tt}
    Simulation(model::Tp, T::Tt, Δt::Ts, S_high::Tsh, S_low::Tsl, t::AbstractVector{Ts}) where {Tp,Tt,Ts,Tsh,Tsl} = Δt > T ? error("The sampling period Δt must be less than the duration T") : new{Tp,Tt,Ts,Tsh,Tsl}(model, T, Δt, S_high, S_low, t)
end

function _regular_sampling(model::PowerSpectralDensity, T::Tt, Δt::Ts, S_high::Tsh, S_low::Tsl) where {Tt<:Real,Ts<:Real,Tsh<:Real,Tsl<:Real}
    t = range(start=0, stop=T, step=Δt)
    return Simulation(model, T, Δt, S_high, S_low, t)
end

Simulation(model::Tp, T::Tt, Δt::Ts) where {Tp<:PowerSpectralDensity,Tt<:Real,Ts<:Real} = _regular_sampling(model, T, Δt, 10.0, 10.0)
Simulation(model::Tp, T::Tt, Δt::Ts, S_high::Tsh, S_low::Tsl) where {Tp<:PowerSpectralDensity,Tt<:Real,Ts<:Real,Tsh<:Real,Tsl<:Real} = _regular_sampling(model, T, Δt, S_high, S_low)


function _arbitrary_sampling(model::PowerSpectralDensity, t::AbstractVector{Tt}, S_high::Tsh, S_low::Tsl) where {Tt<:Real,Tsh<:Real,Tsl<:Real}
    if sort(t) != t
        error("The time vector must be sorted")
    end
    T = t[end] - t[1]
    Δt = minimum(diff(t))
    return Simulation(model, T, Δt, S_high, S_low, t)
end

Simulation(model::Tp, t::AbstractVector{Tt}) where {Tp<:PowerSpectralDensity,Tt<:Real} = _arbitrary_sampling(model, t, 10.0, 10.0)
Simulation(model::Tp, t::AbstractVector{Tt}, S_high::Tsh, S_low::Tsl) where {Tp<:PowerSpectralDensity,Tt<:Real,Tsh<:Real,Tsl<:Real} = _arbitrary_sampling(model, t, S_high, S_low)

"""
    timmer_koenig(psd, rng)

Generate a time series with a given power spectral density (PSD) using the Timmer & Koenig method.

# Arguments
- `psd::Array{Float64, 1}`: PSD of the time series.
- `rng::MersenneTwister`: Random number generator.

# Returns
- `x::Array{Float64, 1}`: Time series with the given PSD.
"""
function timmer_koenig(psd, rng::Random.AbstractRNG)
    N = length(psd)
    Num = randn(rng, (N, 2))
    Re, Im = Num[:, 1], Num[:, 2]

    Im[end] = 0.0
    Rand_psd = sqrt.(psd / 2) .* (Re + im * Im)
    insert!(Rand_psd, 1, 0.0)# N+1 frequencies including 0 and Nyquist
    x = irfft(Rand_psd, 2 * (N + 1) - 2) # 2*(N+1)-2 is the length of the time series
    return x[1:N+1] * N
end

"""
    timmer_koenig_alt(psd, rng)

Generate a time series with a given power spectral density (PSD) using the Timmer & Koenig method.

# Arguments
- `psd::Array{Float64, 1}`: PSD of the time series.
- `rng::MersenneTwister`: Random number generator.

# Returns
- `x::Array{Float64, 1}`: Time series with the given PSD.
"""
function timmer_koenig_alt(psd, rng::Random.AbstractRNG)
    N = length(psd)
    Χ²₂ = rand(rng, Exponential(2), N - 1)
    Χ²₁ = rand(rng, Chisq(1), 1)
    θ = rand(rng, N)
    θ[end] = 0.0
    Rand_psd = sqrt.(psd / 2 .* vcat(Χ²₂, Χ²₁)) .* exp.(2π * im * θ)

    insert!(Rand_psd, 1, 0.0)# N+1 frequencies including 0 and Nyquist
    x = irfft(Rand_psd, 2 * (N + 1) - 2) # 2*(N+1)-2 is the length of the time series
    return x[1:N+1] * N
end

"""
    sample(rng, sim, n)

Generate a time series with a given power spectral density (PSD) using the Timmer & Koenig method.

# Arguments
- `rng::MersenneTwister`: Random number generator.
- `sim::Simulation`: The simulation 
- `n::Int`: Number of time series to generate.

"""
function Distributions.sample(rng::Random.AbstractRNG, sim::Simulation, n::Int=1, input_mean=0; Fvar=nothing, alt::Bool=false, poisson=false, exponentiate=false, error_size=0.02)
    Δf = 1 / sim.T / sim.S_low
    fₘ = 1 / sim.Δt / 2 * sim.S_high
    f = range(start=Δf, step=Δf, stop=fₘ)

    psd = sim.model(f)
    # get the "true" time series
    t = range(start=0, step=1 / (2fₘ), stop=1 / (2Δf))
    if alt
        x = [timmer_koenig_alt(psd, rng) * sqrt(Δf) for i in 1:n]
    else
        x = [timmer_koenig(psd, rng) * sqrt(Δf) for i in 1:n]
    end
    x = hcat(x...)

    # add the mean
    xm = mean(x, dims=1)
    xstd = std(x, dims=1)


    # change the time series to the desired time vector
    indexes = searchsortedlast.(Ref(t), sim.t)
    xₛ = x[indexes, :]


    if !isnothing(Fvar)
        xₛ = ((xₛ .- xm) ./ xstd * Fvar  .+ 1 ) .* input_mean
    else
        xₛ = (xₛ .- xm) .+ input_mean
    end
    # randomise the time series
    if poisson
        if any(xₛ .<= 0)
            println("Warning: Poisson noise is only valid for positive values. Setting to 0.")
            xₛ[xₛ.<=0] .= 0
        end
        x = rand.(rng, Poisson.(xₛ .* sim.Δt)) ./ sim.Δt
        σₓ = sqrt.(x .* sim.Δt) ./ sim.Δt

    else
        if exponentiate
            xₛ = exp.(xₛ)
        end
        σₓ = sqrt.(abs.(xₛ)) * error_size .* abs.(randn(rng, size(xₛ)))
        x = xₛ + σₓ .* randn(rng, size(xₛ))
    end

    return simu.t, x, σₓ
end

Distributions.sample(sim::Simulation, n::Int=1, alt::Bool=false) = Distributions.sample(Random.GLOBAL_RNG, sim, n, alt)
