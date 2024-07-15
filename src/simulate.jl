@doc raw"""
	Simulation(model::PowerSpectralDensity, T::Real, Δt::Real, S_high::Real, S_low::Real, t::AbstractVector{Real})

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
struct Simulation{Tp <: PowerSpectralDensity, Tt <: Real, Ts <: Real, Tsh <: Real, Tsl <: Real}
	model::Tp
	T::Tt
	Δt::Ts
	S_high::Tsh
	S_low::Tsl
	t::AbstractVector{Tt}
	Simulation(model::Tp, T::Tt, Δt::Ts, S_high::Tsh, S_low::Tsl, t::AbstractVector{Ts}) where {Tp, Tt, Ts, Tsh, Tsl} = Δt > T ? error("The sampling period Δt must be less than the duration T") : new{Tp, Tt, Ts, Tsh, Tsl}(model, T, Δt, S_high, S_low, t)
end

function _regular_sampling(model::PowerSpectralDensity, T::Tt, Δt::Ts, S_high::Tsh, S_low::Tsl) where {Tt <: Real, Ts <: Real, Tsh <: Real, Tsl <: Real}
	t = range(start = 0, stop = T, step = Δt)
	return Simulation(model, T, Δt, S_high, S_low, t)
end

Simulation(model::Tp, T::Tt, Δt::Ts) where {Tp <: PowerSpectralDensity, Tt <: Real, Ts <: Real} = _regular_sampling(model, T, Δt, 10.0, 10.0)
Simulation(model::Tp, T::Tt, Δt::Ts, S_high::Tsh, S_low::Tsl) where {Tp <: PowerSpectralDensity, Tt <: Real, Ts <: Real, Tsh <: Real, Tsl <: Real} = _regular_sampling(model, T, Δt, S_high, S_low)


function _arbitrary_sampling(model::PowerSpectralDensity, t::AbstractVector{Tt}, S_high::Tsh, S_low::Tsl) where {Tt <: Real, Tsh <: Real, Tsl <: Real}
	if sort(t) != t
		error("The time vector must be sorted")
	end
	T = t[end] - t[1]
	Δt = minimum(diff(t))
	return Simulation(model, T, Δt, S_high, S_low, t)
end

Simulation(model::Tp, t::AbstractVector{Tt}) where {Tp <: PowerSpectralDensity, Tt <: Real} = _arbitrary_sampling(model, t, 10.0, 10.0)
Simulation(model::Tp, t::AbstractVector{Tt}, S_high::Tsh, S_low::Tsl) where {Tp <: PowerSpectralDensity, Tt <: Real, Tsh <: Real, Tsl <: Real} = _arbitrary_sampling(model, t, S_high, S_low)

@doc raw"""
	timmer_koenig(psd, rng::random.AbstractRNG, α=1.0)

Generate a time series with a given power spectral density (PSD) using the [1995A&A...300..707T](@cite) method

1. Given N values of the power spectral density (PSD) 𝓟
2. Draw 2N values from a standard normal distribution. 
3. The amplitude of the randomised periodogram is given by A = N + i M
3. The randomised periodogram is given by 𝓟_rand = √(𝓟 / 2) * A
5. The first value of the randomised periodogram inserted and set to 0
6. The time series is obtained by taking the inverse Fourier transform of the randomised periodogram, with the function FFTW.irfft. We use the length of the time series as 2*(N+1)-1.

# Arguments
- `𝓟::Array{Float64, 1}`: PSD associated with the process
- `rng::Random.AbstractRNG`: Random number generator.
- `α::Real`: Multiplicative factor for the randomised periodogram. Default is 1.0. This can be a complex vector to add a phase to the time series.

# Returns
- `x::Array{Float64, 1}`: Time series with the given PSD.
"""
function timmer_koenig(𝓟, rng::Random.AbstractRNG)
	N = length(𝓟)
	Num = randn(rng, (N, 2))
	Re, Im = Num[:, 1], Num[:, 2]

	Im[end] = 0.0
	Rand_psd = .√(𝓟 / 2) .* (Re + im * Im)
	insert!(Rand_psd, 1, 0.0)# N+1 frequencies including 0 and Nyquist
	x = irfft(Rand_psd, 2 * (N + 1) - 1) # 2*(N+1)-2 is the length of the time series
	return x
end

@doc raw"""
	timmer_koenig_alt(𝓟, rng::Random.AbstractRNG, α = 1.0)

Generate a time series with a given power spectral density (PSD) using the [1995A&A...300..707T](@cite) method with an alternative parametrisation

1. Given N values of the power spectral density (PSD) 𝓟
2. Draw N-1 values from a χ²(2)  distribution and 1 value from χ²₁(1), this is A the amplitude of the randomised periodogram
3. Draw N values from a uniform distribution between 0 and 1, and set the last value to 0, this is θ the phase of the randomised periodogram
4. The randomised periodogram is given by 𝓟_rand = √(𝓟 / 2 * A) * exp(2πiθ)
5. The first value of the randomised periodogram inserted and set to 0
6. The time series is obtained by taking the inverse Fourier transform of the randomised periodogram, with the function FFTW.irfft. We use the length of the time series as 2*(N+1)-1.

# Arguments
- `𝓟::Array{Float64, 1}`: PSD associated with the process
- `rng::Random.AbstractRNG`: Random number generator.
- `α::Real`: Multiplicative factor for the randomised periodogram. Default is 1.0. This can be a complex vector to add a phase to the time series.

# Returns
- `x::Array{Float64, 1}`: Time series with the given PSD.
"""
function timmer_koenig_alt(𝓟, rng::Random.AbstractRNG, α = 1.0)
	N = length(𝓟)
	Χ²₂ = rand(rng, Exponential(2), N - 1)
	Χ²₁ = rand(rng, Chisq(1), 1)
	θ = rand(rng, N)
	θ[end] = 0.0
	Rand_psd = .√(𝓟 / 2 .* vcat(Χ²₂, Χ²₁)) .* exp.(2π * im * θ) .* α

	insert!(Rand_psd, 1, 0.0)# N+1 frequencies including 0 and Nyquist
	x = irfft(Rand_psd, 2 * (N + 1) - 1) # 2*(N+1)-2 is the length of the time series
	return x
end

@doc raw"""
Split a long time series into shorter time series.

Break the time series into `n_slices` shorter time series. The short time series are of equal length.

# Arguments
- `t`: The time indexes of the long time series.
- `ts`: The values of the long time series.
- `n_slices`: The number of slices to break the time series into.

# Returns
- A tuple of two lists: the first containing the time indexes of the shorter time series, and the second containing the values of the shorter time series.
"""
function split_longtimeseries(t, ts, n_slices::Int)
	t_slices = []
	ts_slices = []
	size_slice = div(length(t), n_slices)
	for i in 1:n_slices
		start_index = (i - 1) * size_slice + 1
		end_index = i * size_slice
		t_slice = collect(t[start_index:end_index])
		ts_slice = ts[start_index:end_index]

		push!(t_slices, t_slice)
		push!(ts_slices, ts_slice)
	end
	return t_slices, ts_slices
end

@doc raw"""
	sample(rng, sim, n=1, input_mean=0; split_long=false, randomise_values=true, Fvar=nothing, alt=false, poisson=false, exponentiate=false, error_size=0.02)

Generate a time series with a given power spectral density (PSD) using the [1995A&A...300..707T](@cite) method.

# Arguments
- `rng::MersenneTwister`: Random number generator.
- `sim::Simulation`: The simulation 
- `n::Int`: Number of time series to generate. Default is 1.
- `input_mean::Real`: The mean of the time series. Default is 0.
- `split_long::Bool`: If true, the time series is split into shorter time series given by `sim.S_low`. Default is false.
- `randomise_values::Bool`: If true, the values of the time series are randomised. Default is true.
- `Fvar::Real`: The variance of the time series. Default is nothing.
- `alt::Bool`: If true, uses the alternative Timmer & Koenig method. Default is false. 
- `poisson::Bool`: If true, Poisson noise is added to the time series. Default is false.
- `exponentiate::Bool`: If true, the time series is exponentiated. Default is false.
- `error_size::Real`: The size of the error. Default is 0.02.

"""
function Distributions.sample(rng::Random.AbstractRNG, sim::Simulation, n::Int = 1, input_mean = 0; randomise_values=true, split_long = false, Fvar = nothing, alt::Bool = false, poisson = false, exponentiate = false, error_size = 0.02)
	Δf = 1 / sim.T / sim.S_low
	fₘ = 1 / sim.Δt / 2 * sim.S_high
	Δτ = 1 / 2fₘ
	f = range(start = Δf, step = Δf, stop = fₘ)

    if sim.model isa PowerSpectralDensity
        psd = sim.model(f)
        # get the "true" time series
        if alt
            x = [timmer_koenig_alt(psd, rng) * sqrt(2Δf) * length(f) for i in 1:n]
        else
            x = [timmer_koenig(psd, rng) * sqrt(2Δf) * length(f) for i in 1:n]
        end
        x = hcat(x...)
        t = range(0, step = Δτ, length = size(x, 1))
    end

	n_slices = round(Int,sim.S_low)

    # split the long time series
	if split_long
		if n_slices < 5
			println("Warning: The number of slices is less than 5. Setting S_low to a high value is required.")
			println("Not splitting the time series.")
			indexes = searchsortedlast.(Ref(t), sim.t)
			xₛ = x[indexes, :]
		else
			xₛ = []
			if n == 1 # for a single long time series
				t_long, x_long = split_longtimeseries(t, x, n_slices)

				for i in 1:n_slices
					tcurr = collect(t_long[i]) .- t_long[i][1]
					indexes = searchsortedlast.(Ref(tcurr), sim.t)
					push!(xₛ, x_long[i][indexes])
				end

			else
				for j in 1:n # for set of long time series
					t_long, x_long = split_longtimeseries(t, reshape(x[:, j], (size(x, 1), 1)), n_slices)

					for i in 1:n_slices
						tcurr = collect(t_long[i]) .- t_long[i][1]
						indexes = searchsortedlast.(Ref(tcurr), sim.t)
						push!(xₛ, x_long[i][indexes])
					end
				end
			end
			xₛ = hcat(xₛ...)
			t = t_long[1]
		end
	else
        # change the time series to the desired time vector
		indexes = searchsortedlast.(Ref(t), sim.t)
		xₛ = x[indexes, :]
	end
    times = t[indexes]

    # return the time series without randomising the values
    if !randomise_values
        return times, xₛ, zeros(size(xₛ))
    end

	# add the mean
	xm = mean(xₛ, dims = 1)
	xstd = std(xₛ, dims = 1)

	if !isnothing(Fvar)
		xₛ = ((xₛ .- xm) ./ xstd * Fvar .+ 1) .* input_mean
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

	return times, x, σₓ
end

Distributions.sample(sim::Simulation, n::Int = 1, alt::Bool = false) = Distributions.sample(Random.GLOBAL_RNG, sim, n, alt)
