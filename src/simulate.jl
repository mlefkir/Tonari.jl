@doc raw"""
	Simulation(model::Model, T::Real, Δt::Real, S_high::Real, S_low::Real, t::AbstractVector{Real})

A struct that contains the information of a simulation of a stochastic process.

# Fields
- `model::Model`: The power spectral density of the process.
- `T::Real`: The duration of the simulated time series. Note that the time stamps are from 0 to `T - Δt`, so the duration is `T - Δt`.
- `Δt::Real`: The sampling period, or minimum time difference between samples.
- `S_high::Real`: The factor by which the maximum frequency is multiplied for the simulation.
- `S_low::Real`: The factor by which the minimum frequency is divided for the simulation.
- `t::AbstractVector{Real}`: The time vector, in this case, the time at which the process is sampled. It is assumed that the time vector is sorted.

# Constructors
- `Simulation(model::Model, T::Real, Δt::Real, S_high::Real, S_low::Real, t::AbstractVector{Real})`: Constructs a simulation with regular sampling.
- `Simulation(model::Model, T::Real, Δt::Real)`: Constructs a simulation with regular sampling, and sets `S_high` and `S_low` to 10.0.
- `Simulation(model::Model, T::Real, Δt::Real, S_high::Real, S_low::Real)`: Constructs a simulation with the given parameters.
- `Simulation(model::Model, t::AbstractVector{Real}, S_high::Real, S_low::Real)`: Constructs a simulation with sampling pattern given by `t`.
- `Simulation(model::Model, t::AbstractVector{Real})`: Constructs a simulation with sampling pattern given by `t`, and sets `S_high` and `S_low` to 10.0.
"""
struct Simulation{Tp <: Model, Tt <: Real, Ts <: Real, Tsh <: Real, Tsl <: Real}
	model::Tp
	T::Tt
	Δt::Ts
	S_high::Tsh
	S_low::Tsl
	t::AbstractVector{Tt}
	Simulation(model::Tp, T::Tt, Δt::Ts, S_high::Tsh, S_low::Tsl, t::AbstractVector{Ts}) where {Tp, Tt, Ts, Tsh, Tsl} = Δt > T ? error("The sampling period Δt must be less than the duration T") : new{Tp, Tt, Ts, Tsh, Tsl}(model, T, Δt, S_high, S_low, t)
end

function _regular_sampling(model::Model, T::Tt, Δt::Ts, S_high::Tsh, S_low::Tsl) where {Tt <: Real, Ts <: Real, Tsh <: Real, Tsl <: Real}
	t = range(start = 0, stop = T , step = Δt)
	return Simulation(model, T , Δt, S_high, S_low, t)
end

Simulation(model::Tp, T::Tt, Δt::Ts) where {Tp <: Model, Tt <: Real, Ts <: Real} = _regular_sampling(model, T, Δt, 10.0, 10.0)
Simulation(model::Tp, T::Tt, Δt::Ts, S_high::Tsh, S_low::Tsl) where {Tp <: Model, Tt <: Real, Ts <: Real, Tsh <: Real, Tsl <: Real} = _regular_sampling(model, T, Δt, S_high, S_low)


function _arbitrary_sampling(model::Model, t::AbstractVector{Tt}, S_high::Tsh, S_low::Tsl) where {Tt <: Real, Tsh <: Real, Tsl <: Real}
	if !issorted(t)
		error("The time vector must be sorted")
	end
	T = t[end] - t[1]
	Δt = minimum(diff(t))
	return Simulation(model, T, Δt, S_high, S_low, t .- t[1])
end

Simulation(model::Tp, t::AbstractVector{Tt}) where {Tp <: Model, Tt <: Real} = _arbitrary_sampling(model, t, 10.0, 10.0)
Simulation(model::Tp, t::AbstractVector{Tt}, S_high::Tsh, S_low::Tsl) where {Tp <: Model, Tt <: Real, Tsh <: Real, Tsl <: Real} = _arbitrary_sampling(model, t, S_high, S_low)

@doc raw"""
	timmer_koenig(psd, rng::random.AbstractRNG, alternative=false)

Generate a time series with a given power spectral density (PSD) using the [1995A&A...300..707T](@cite) method

1. Given N values of the power spectral density (PSD) 𝓟
2. Draw 2N values from a standard normal distribution. 
3. The amplitude of the randomised periodogram is given by A = N + i M
3. The randomised periodogram is given by 𝓟_rand = √(𝓟 / 2) * A
5. The first value of the randomised periodogram inserted and set to 0
6. The time series is obtained by taking the inverse Fourier transform of the randomised periodogram, with the function FFTW.irfft. We use the length of the time series as 2*(N+1)-1.

The alternative parametrisation is given by

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
function timmer_koenig(𝓟, rng::Random.AbstractRNG; alternative = false)
	N = length(𝓟)
	Rand_psd = get_randomised_psd(𝓟, rng, alternative)

	insert!(Rand_psd, 1, 0.0)# N+1 frequencies including 0 and Nyquist
	x = irfft(Rand_psd, 2 * (N + 1) - 1) # 2*(N+1)-1 is the length of the time series
	return x
end

function timmer_koenig_manual(Rand_psd)
	N = length(Rand_psd)
	insert!(Rand_psd, 1, 0.0)# N+1 frequencies including 0 and Nyquist
	x = irfft(Rand_psd, 2 * (N + 1) - 1) # 2*(N+1)-1 is the length of the time series
	return x
end

function get_randomised_psd(𝓟, rng::Random.AbstractRNG, alternative = false)
	N = length(𝓟)
	if alternative
		Num = randn(rng, (N, 2))
		Re, Im = Num[:, 1], Num[:, 2]
		Im[end] = 0.0
		Rand_psd = .√(𝓟 / 2) .* (Re + im * Im)
	else
		Χ²₂ = rand(rng, Exponential(2), N - 1)
		Χ²₁ = rand(rng, Chisq(1), 1)
		θ = rand(rng, N)
		θ[end] = 0.0
		Rand_psd = .√(𝓟 / 2 .* vcat(Χ²₂, Χ²₁)) .* exp.(2π * im * θ)
	end
	return Rand_psd
end

@doc raw"""
	split_longtimeseries(t, ts, n_slices, t_end)

Split a long time series into shorter time series.
Break the time series into `n_slices` shorter time series. The short time series are of equal length.

# ArgumentsT
- `t`: The time indexes of the long time series.
- `ts`: The values of the long time series.
- `n_slices`: The number of slices to break the time series into.
- `t_end`

# Returns
- A tuple of two lists: the first containing the time indexes of the shorter time series, and the second containing the values of the shorter time series.
"""
function split_longtimeseries(t, ts, n_slices::Int, t_end)
	t_slices = []
	ts_slices = []
	end_long = sum(t .<= t_end)

	indexes = [i*end_long+1:(end_long)*(i+1) for i in 0:n_slices-1]
	for i in 1:n_slices
		t_slice = t[indexes[i]] .- t[indexes[i]][1]
		ts_slice = ts[indexes[i]]

		push!(t_slices, t_slice)
		push!(ts_slices, ts_slice)
	end
	return t_slices, ts_slices
end


@doc raw"""
	findnearest(a, b)

Find the nearest value in `b` to each value in `a`.
This assumes that a and b are sorted.
"""
function findnearest(a, b)
	@assert issorted(a) && issorted(b)
	indexes = []
	j = 1
	for i in eachindex(a)
		while !(b[j] ≈ a[i]) && j < length(b)
			j += 1
		end
		push!(indexes, j)
	end
	return indexes
end

@doc raw"""
	sample_timeseries(t, y, M)

Extract a random subset of points from the time series.
"""
function sample_timeseries(t, t_desired)
	findnearest(t_desired, t)
end

@doc raw"""
	sample_split_timeseries(x, t, t_desired, n_sim, n, n_slices, split_long)

Split the time series into shorter `.
"""
function sample_split_timeseries(x, t, t_desired, n_sim, n, n_slices, split_long)

	t_end = t_desired[end]
	n_bands = 1
	xₛ = []
	if x isa Vector{Matrix{Float64}}
		# there are multiple bands of time series
		n_bands = length(x)

	end
	if split_long
		if n_slices < 5
			@warn "The number of slices is less than 5. Setting S_low to a high value is required.\nNot splitting the time series."
			indexes = sample_timeseries(t, t_desired)
			for j in 1:n_bands
				if n_bands == 1
					push!(xₛ, x[indexes, :])
				else
					push!(xₛ, x[j][indexes, :])
				end
			end

		else
			for k in 1:n_bands
				push!(xₛ, [])
				l = 0

				for j in 1:n_sim # for set of long time series
					if n_bands == 1
						t_long, x_long = split_longtimeseries(t, x, n_slices, t_end)
					else
						t_long, x_long = split_longtimeseries(t, reshape(x[k][:, j], (size(x[k], 1), 1)), n_slices, t_end)
					end
					for i in 1:n_slices
						tcurr = collect(t_long[i]) .- t_long[i][1]
						indexes = sample_timeseries(tcurr, t_desired)
						push!(xₛ[k], x_long[i][indexes])
						l += 1
						if l == n # if we have the desired number of time series
							break
						end
					end
					times = t_long[1][indexes]
				end
			end

			xₛ_full = [hcat(xₛ[j]...) for j in 1:n_bands]
		end
	else
		# get a subset of the time series
		indexes_t = t .<= (t_desired[end] + √eps(Float64)) # just in case the last value is not included due to floating point errors
		t = t[indexes_t]
		# x = x[:][indexes_t, :]
		# change the time series to the desired time vector
		indexes = sample_timeseries(t, t_desired)
		for j in 1:n_bands
			if n_bands == 1
				push!(xₛ, x[indexes_t, :][indexes, :])
			else
				push!(xₛ, x[j][indexes_t, :][indexes, :])
			end
		end
		xₛ_full = xₛ
		times = t[indexes]
	end
	if n_bands == 1
		xₛ_full = xₛ_full[1]
	end
	@assert t_desired[1:end-1] ≈ times[1:end-1] "The sampled times are not approximatively equal to the desired times"
	if t_desired[end]≈times[end]
		return times, xₛ_full
	else
		@warn "Removing the last point as it is not at the right time stamp"
		#remove the last point
		y = []
		for j in 1:n_bands
			if n_bands == 1
				y = xₛ_full[1:end-1,:]
				# push!(y,xₛ_full[1:end-1,:])
			else
				push!(y,xₛ_full[j][1:end-1,:])
			end
		end
		return times[1:end-1],y
	end
end

@doc raw"""
	draw_errorbars(rng, x, Δt, poisson=false, Fvar=nothing, error_size=0.05)

Draw errorbars for the time series.

# Arguments
- `rng::MersenneTwister`: Random number generator.
- `x::Array{Float64, 1}`: The values of the time series.
- `Δt::Real`: The sampling period.
- `poisson::Bool`: If true, Poisson noise is added to the time series. Default is false.
- `Fvar::Real`: The variance of the time series. Default is nothing.
- `error_size::Real`: The size of the error. Default is 0.05.
"""
function draw_errorbars(rng, x, Δt, poisson = false, Fvar = nothing, error_size = 0.05)

	if poisson
		σₓ = sqrt.(x .* Δt) ./ Δt
	else
		σₓ = sqrt.(abs.(x)) * error_size .* abs.(randn(rng, size(x)))
	end
	return σₓ
end

@doc raw"""
	randomise_fluxes(rng, xₛ, Δt, input_mean = 0.0; σ = nothing, Fvar = nothing, poisson = false, exponentiate = false, error_size = 0.05)

Randomise the fluxes of the time series.

# Arguments
- `rng::MersenneTwister`: Random number generator.
- `xₛ::Array{Float64, 1}`: The values of the time series.
- `Δt::Real`: The sampling period.
- `input_mean::Real`: The mean of the time series. Default is 0.
- `σ::Array{Float64, 1}`: The errorbars of the time series. Default is nothing.
- `Fvar::Real`: The variance of the time series. Default is nothing.
- `poisson::Bool`: If true, Poisson noise is added to the time series. Default is false.
- `exponentiate::Bool`: If true, the time series is exponentiated. Default is false.
"""
function randomise_fluxes(rng, xₛ, Δt, input_mean = 0.0; σ = nothing, Fvar = nothing, poisson = false, exponentiate = false, error_size = 0.05)
	# add the mean
	if !isnothing(σ) && poisson
		@warn "Errorbars given but Poisson distribution chosen, discarding the errorbars"
	end
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
			@warn "Poisson noise is only valid for positive values. Setting to 0."
			xₛ[xₛ.<=0] .= 0
		end
		x = rand.(rng, Poisson.(xₛ .* Δt)) ./ Δt
		σₓ = draw_errorbars(rng, x, Δt, true, Fvar, error_size)

	else
		if exponentiate
			xₛ = exp.(xₛ)
		end
		if isnothing(σ)
			σₓ = draw_errorbars(rng, xₛ, Δt, false, Fvar, error_size)
		else
			@assert all(σ .>= 0.0) "The errorbars are not positive"
			σₓ = σ
		end
		x = xₛ + σₓ .* randn(rng, size(xₛ))
		if exponentiate && any(x .<= 0)
			@warn "Exponentiated time series has negative values. Setting to 0."
			x[x.<=0] .= 0
		end
	end
	return x, σₓ
end

@doc raw"""
	sample(rng, sim, n=1, input_mean=0; split_long=false, randomise_values=true, Fvar=nothing, alt=false, poisson=false, exponentiate=false, error_size=0.02)

Generate a time series with a given power spectral density (PSD) using the [1995A&A...300..707T](@cite) method.

# Arguments
- `rng::MersenneTwister`: Random number generator.
- `sim::Simulation`: The simulation 
- `n::Int`: Number of time series to generate. Default is 1.
- `input_mean::Real`: The mean of the time series. Default is 0.
- `split_long::Bool`: If true, the time series is split into shorter time series given by `sim.S_low`-1. Default is true.
- `randomise_values::Bool`: If true, the values of the time series are randomised. Default is true.
- `Fvar::Real`: The variance of the time series. Default is nothing.
- `alt::Bool`: If true, uses the alternative Timmer & Koenig method. Default is false. 
- `poisson::Bool`: If true, Poisson noise is added to the time series. Default is false.
- `exponentiate::Bool`: If true, the time series is exponentiated. Default is false.
- `error_size::Real`: The size of the error. Default is 0.05.

"""
function Distributions.sample(rng::Random.AbstractRNG, sim::Simulation, n::Int = 1, input_mean = 0; σₓ = nothing, randomise_values = true, split_long = true, Fvar = nothing, alt::Bool = false, poisson = false, exponentiate = false, error_size = 0.05)
	@assert n > 0 "The number of simulations n must be a postive integer!, n=$n is not accepted"
	if !isnothing(Fvar)
		@assert Fvar > 0 "The variance Fvar must be a positive number!"
	end
	if !isnothing(σₓ)
		@assert all(σₓ .> 0) "The errorbars given must be positve!"
	end
	if poisson && exponentiate
		@warn "Poisson noise is not valid for exponentiated time series. Setting poisson to false."
		poisson = false
	end

	Δf = 1 / sim.T / sim.S_low
	fₘ = 1 / sim.Δt / 2 * sim.S_high
	Δτ = 1 / 2fₘ
	f = range(start = Δf, step = Δf, stop = fₘ + Δf)

	n_slices = round(Int, sim.S_low) - 1 # sometimes the long time series cannot be split into S_low slices so we need to adjust this value
	# adjust n if the time series is split
	if split_long
		n_sim = ceil(Int, n / n_slices)
	else
		n_sim = n
	end


	if sim.model isa PowerSpectralDensity
		psd = sim.model(f)
		# get the "true" time series
		x = [timmer_koenig(psd, rng, alternative = alt) * sqrt(2Δf) * length(f) for i in 1:n_sim]
		x = hcat(x...)

		t = 0:Δτ:(size(x, 1)-1)*Δτ

	elseif sim.model isa CrossSpectralDensity

		Δφ = exp.(im * 2π .* f .* calculate(sim.model.Δφ, f))

		psd = sim.model.𝓟₁(f)
		# get the randomised periodogram
		X = [get_randomised_psd(psd, rng, alt) * sqrt(2Δf) * length(f) for i in 1:n_sim]
		X = hcat(X...)

		if sim.model.𝓟₁ == sim.model.𝓟₂
			x₁ = [timmer_koenig_manual(X[:, i]) for i in 1:n_sim]
			x₂ = [timmer_koenig_manual(X[:, i] .* Δφ) for i in 1:n_sim]
		else
			psd₁ = sim.model.𝓟₁(f)
			psd₂ = sim.model.𝓟₂(f)
			ratio = sqrt.(psd₂ ./ psd₁)
			x₁ = [timmer_koenig_manual(X[:, i]) for i in 1:n_sim]
			x₂ = [timmer_koenig_manual(X[:, i] .* ratio .* Δφ) for i in 1:n_sim]
		end

		x₁ = hcat(x₁...)
		x₂ = hcat(x₂...)
		x = [x₁, x₂]
		t = range(0, step = Δτ, length = size(x₁, 1))
	else
		error("The model must be a PowerSpectralDensity or CrossSpectralDensity")
	end

	# split the long time series and resample at the desired time stamps
	times, xₛ = sample_split_timeseries(x, t, sim.t[1:end-1], n_sim, n, n_slices, split_long)

	# return the time series without randomising the values
	if !randomise_values
		return times, xₛ, zeros(size(xₛ))
	end

	# randomise the fluxes
	if sim.model isa PowerSpectralDensity
		x, σₓ = randomise_fluxes(rng, xₛ, sim.Δt, input_mean, σ = σₓ, Fvar = Fvar, poisson = poisson, exponentiate = exponentiate, error_size = error_size)
	elseif sim.model isa CrossSpectralDensity
		x, σₓ = [], []
		for i in 1:2
			x_band, σₓ_band = randomise_fluxes(rng, xₛ[i], sim.Δt, input_mean, Fvar = Fvar, poisson = poisson, exponentiate = exponentiate, error_size = error_size)
			push!(x, x_band)
			push!(σₓ, σₓ_band)
		end
	else
		error("The model must be a PowerSpectralDensity or CrossSpectralDensity")
	end

	return times, x, σₓ
end

Distributions.sample(sim::Simulation, n::Int = 1, alt::Bool = false) = Distributions.sample(Random.GLOBAL_RNG, sim, n, alt)
