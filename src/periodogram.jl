"""
	end_matching(y,t)

End-match the data y with a straight line and subtract it from the data.
"""
function end_matching(y, t)
	a = (y[end] - y[1]) / (t[end] - t[1])
	b = y[1]
	y_detrended = y - (a * t .+ b)
	return y_detrended
end

@doc raw"""
	periodogram(y, t, normalisation = "default"; apply_end_matching = false, subtract_mean = true)

Compute the periodogram of the data y with time stamps t using the fast Fourier transform (FFT).

The periodogram is computed as the squared magnitude of the Fourier transform of the data. In practrice we 
use the real-valued fast Fourier transform (rfft) to compute the periodogram.

# Arguments
- `y::Array{Float64, 1}`: Time series data.
- `t::Array{Float64, 1}`: Time stamps.
- `normalisation::Float64`: Normalisation factor. Default is 2Δt / length(t).
- `apply_end_matching::Bool`: Apply end-matching to the data. Default is false.
- `subtract_mean::Bool`: Subtract the mean from the data. Default is true.
"""
function periodogram(t, y, normalisation = "default"; apply_end_matching = false, subtract_mean = true)

	@assert length(t) == size(y, 1) "The length of the time stamps must be the same as the length of the time series."

	average_periodogram = false
	if ndims(y) == 2
		if size(y, 2) > 1
			println("Multiple time series detected. Computing the average periodogram.")
			average_periodogram = true
		end
	end

	# define the frequency range
	Δt = t[2] - t[1]
	T = t[end] - t[1]
	fmin, fmax = 1 / T, 1 / (2Δt)
	# f = range(fmin, step=fmin, fmax)
	f = rfftfreq(size(y, 1), 1 / Δt)

	if normalisation == "default"
		normalisation = 2Δt / length(t)
	else
		error("Normalisation $normalisation not recognised.")
	end

	# detrend the data
	if apply_end_matching
		if average_periodogram
			x = [end_matching(y[:, i], t) for i in range(1, size(y, 2))]
			x = hcat(x...)
		else
			x = end_matching(y, t)
		end
	else
		x = y
	end

	# subtract the mean from the data
	if subtract_mean
		if average_periodogram
			P = [abs.(rfft(x[:, i] .- mean(x[:, i]))) .^ 2 for i in range(1, size(x, 2))]
			P = sum(P, dims = 1)[1] / size(x, 2)
		else
			P = abs.(rfft(x .- mean(x))) .^ 2
		end
	else
		if average_periodogram
			P = [abs.(rfft(x[:, i])) .^ 2 for i in range(1, size(x, 2))]
			P = sum(P, dims = 1)[1] / size(x, 2)
		else
			P = abs.(rfft(x)) .^ 2
		end
	end

	# check that the length of the periodogram is the same as the frequency range
	if length(f) != length(P)
		error("The length of the periodogram $(length(P)) is not the same as the frequency range $(length(f)).")
	end
	return f[2:end], P[2:end] * normalisation
end

@doc raw"""
	cross_periodogram(y₁, y₂, t, normalisation = "default"; apply_end_matching = false, subtract_mean = true)

	

"""
function cross_periodogram(t, y₁, y₂, y₁_err = nothing, y₂_err = nothing, normalisation = "default"; compute_coherence = true, apply_end_matching = false, subtract_mean = true)

	# checks on the input data
	@assert size(y₁) == size(y₂) "The size of the two time series must be the same."
	@assert length(t) == size(y₁, 1) "The length of the time stamps must be the same as the length of the time series."
	if !isnothing(y₁_err)
		@assert size(y₁) == size(y₁_err) "The size of the time series and the error must be the same."
	end
	if !isnothing(y₂_err)
		@assert size(y₂) == size(y₂_err) "The size of the time series and the error must be the same."
	end

	if ndims(y₁) == 2
		if size(y₁, 2) == 1
			if compute_coherence
				error("Cannot compute coherence or lags without several segments of the full time series.")
			end
		else
			println("Multiple time series detected. Computing the average cross-periodogram.")
			average_periodogram = true
			n_segments = size(y₁, 2)
			if n_segments < 20
				println("Warning: The number of segments is less than 20. The coherence and lags may not be reliable.")
			end
		end
	end

	# define the frequency range
	Δt = t[2] - t[1]
	T = t[end] - t[1]
	fmin, fmax = 1 / T, 1 / (2Δt)
	f = rfftfreq(size(y₁, 1), 1 / Δt)[2:end] # remove the zero frequency
	f2 = range(fmin, step=fmin, fmax)

	if normalisation == "default"
		normalisation = 2Δt / length(t)
	else
		error("Normalisation $normalisation not recognised.")
	end

	# detrend the data
	if apply_end_matching
		if average_periodogram
			x₁ = [end_matching(y₁[:, i], t) for i in range(1, size(y₁, 2))]
			x₂ = [end_matching(y₂[:, i], t) for i in range(1, size(y₂, 2))]
			x₁, x₂ = hcat(x₁...), hcat(x₂...)
		else
			x₁ = end_matching(y₁, t)
			x₂ = end_matching(y₂, t)
		end
	else
		x₁, x₂ = y₁, y₂
	end

	# compute the DFTs of the two time series
	# subtract the mean from the data
	if subtract_mean
		if average_periodogram
			X₁ = [rfft(x₁[:, i] .- mean(x₁[:, i])) for i in range(1, size(x₁, 2))]
			X₂ = [rfft(x₂[:, i] .- mean(x₂[:, i])) for i in range(1, size(x₂, 2))]
			X₁, X₂ = hcat(X₁...), hcat(X₂...)
			X₁, X₂ = X₁[2:end, :], X₂[2:end, :]
		else
			X₁ = rfft(x₁ .- mean(x₁))
			x₂ = rfft(x₂ .- mean(x₂))
			X₁, X₂ = X₁[2:end], X₂[2:end]
		end
	else
		if average_periodogram
			X₁ = [rfft(x₁[:, i]) for i in range(1, size(x₁, 2))]
			X₂ = [rfft(x₂[:, i]) for i in range(1, size(x₂, 2))]
			X₁, X₂ = hcat(X₁...), hcat(X₂...)
			X₁, X₂ = X₁[2:end, :], X₂[2:end, :]
		else
			X₁ = rfft(x₁)
			X₂ = rfft(x₂)
			X₁, X₂ = X₁[2:end], X₂[2:end]
		end
	end

	# compute the cross-periodogram
	C = conj.(X₁) .* X₂
	C̄ = mean(C, dims = 2)
	# compute individual periodograms
	P₁ = abs.(X₁) .^ 2
	P₂ = abs.(X₂) .^ 2
	P̄₁ = mean(P₁, dims = 2)
	P̄₂ = mean(P₂, dims = 2)

	if compute_coherence
		# compute the raw coherence
		γ² = abs.(C̄) .^ 2 ./ (P̄₁ .* P̄₂)
		γ²_err = √2 .* .√abs.(γ²) .* (1.0 .- γ²) ./ √n_segments # Eq. 9.81 in Bendat and Piersol (2010)

		# compute the phase lags
		Δφ = angle.(C̄)
		Δφ_err = (1 .- γ²) ./ .√(2 .* abs.(γ²) .* n_segments)

		# compute the time lags
		Δτ = Δφ ./ (2π .* f)
		Δτ_err = Δφ_err ./ (2π .* f)

		if !isnothing(y₁_err) && !isnothing(y₂_err)
			# compute corrected coherence

			# noise powers
			N₁ = length(t) * mean(y₁_err .^ 2)
			N₂ = length(t) * mean(y₂_err .^ 2)

			# intrinsic powers
			S₁ = P̄₁ .- N₁
			S₂ = P̄₂ .- N₂

			# 
			n² = (S₁ .* N₂ .+ N₁ .* S₂ .+ N₁ .* N₂) ./ n_segments
			γ²_corrected = (abs.(C̄) .^ 2 .- n²) ./ (S₁ .* S₂)

			a = (2 * n_segments .* n² .^ 2) ./ ((abs.(C̄) .^ 2 .- n²) .^ 2)
			b = N₁^2 ./ S₁ .^ 2 .+ N₂^2 ./ S₂ .^ 2
			c = n_segments .* γ²_err ./ γ² .^ 2
			γ²_corrected_err = 1 / √n_segments .* sqrt.(a .+ b .+ c)
			return f,f2, γ², γ²_corrected, Δφ, γ²_err, γ²_corrected_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, N₁, N₂, n²

		end
		return f, γ², Δφ, γ²_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, C
	end
	return f, C̄, P̄₁, P̄₂, C
end
