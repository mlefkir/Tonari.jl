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
	periodogram(t,y, normalisation = "default"; apply_end_matching = false, subtract_mean = true)

Compute the periodogram of the data y with time stamps t using the fast Fourier transform (FFT).

The periodogram is computed as the squared magnitude of the Fourier transform of the data. In practrice we 
use the real-valued fast Fourier transform (rfft) to compute the periodogram.

# Arguments
- `t::Array{Float64, 1}`: Time stamps.
- `y::Array{Float64, 1}`: Time series data.
- `normalisation::Float64`: Normalisation factor. Default is 2Δt / length(t).
- `apply_end_matching::Bool`: Apply end-matching to the data. Default is false.
- `subtract_mean::Bool`: Subtract the mean from the data. Default is true.
"""
function periodogram(t, y, normalisation = "default"; apply_end_matching = false, subtract_mean = true)

	@assert length(t) == size(y, 1) "The length of the time stamps must be the same as the length of the time series."

	dt = unique(diff(t))
	@assert length(dt) == 1 || all(isapprox.(dt, dt[1])) "The time stamps are not equally spaced! Check your time array!"

	average_periodogram = false
	if ndims(y) == 2
		if size(y, 2) > 1
			@info "Multiple time series detected. Computing the average periodogram."
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
	cross_periodogram(t, y₁, y₂, y₁_err = nothing, y₂_err = nothing; compute_coherence = true, apply_end_matching = false, subtract_mean = true)

Compute the cross-periodogram between two time series y₁ and y₂ with time stamps t using the fast Fourier transform (FFT).

# Arguments
- `t::Array{Float64, 1}`: Time stamps.
- `y₁::Array{Float64, 1}`: First time series data. Or a matrix of multiple time series. 
- `y₂::Array{Float64, 1}`: Second time series data. Or a matrix of multiple time series.
- `y₁_err::Array{Float64, 1}`: Error in the first time series. Or a matrix of errors for multiple time series.
- `y₂_err::Array{Float64, 1}`: Error in the second time series Or a matrix of errors for multiple time series.
- `compute_coherence::Bool`: Compute the coherence between the two time series. Default is true.
- `apply_end_matching::Bool`: Apply end-matching to the data. Default is false.
- `subtract_mean::Bool`: Subtract the mean from the data. Default is true.

If compute_coherence  is true, the function returns the coherence between the two time series. Otherwise, it returns the cross-periodogram and the mean power spectra of the two time series.
see `coherence` function for the optional returns.

# Returns
- `f::Array{Float64, 1}`: Frequency array.
- `C̄::Array{Complex{Float64}, 1}`: Mean cross-periodogram.
- `P̄₁::Array{Float64, 1}`: Mean power spectrum of the first time series.
- `P̄₂::Array{Float64, 1}`: Mean power spectrum of the second time series.
- `C::Array{Complex{Float64}, 1}`: Cross-periodogram.

"""
function cross_periodogram(t, y₁, y₂, y₁_err = nothing, y₂_err = nothing; compute_coherence = true, apply_end_matching = false, subtract_mean = true)

	# checks on the input data
	@assert size(y₁) == size(y₂) "The size of the two time series must be the same."
	@assert length(t) == size(y₁, 1) "The length of the time stamps must be the same as the length of the time series."
	if !isnothing(y₁_err)
		@assert size(y₁) == size(y₁_err) "The size of the time series and the error must be the same."
	end
	if !isnothing(y₂_err)
		@assert size(y₂) == size(y₂_err) "The size of the time series and the error must be the same."
	end
	average_periodogram = false
	if ndims(y₁) == 2
		if size(y₁, 2) == 1
			if compute_coherence
				error("Cannot compute coherence or lags without several segments of the full time series.")
			end
		else
			@info "Multiple time series detected. Computing the average cross-periodogram."
			average_periodogram = true
			n_segments = size(y₁, 2)
			if n_segments < 20
				@warn "Warning: The number of segments is less than 20. The coherence and lags may not be reliable."
			end
		end
	end

	# define the frequency range
	Δt = t[2] - t[1]
	T = t[end] - t[1]
	fmin, fmax = 1 / T, 1 / (2Δt)
	f = rfftfreq(size(y₁, 1), 1 / Δt)[2:end] # remove the zero frequency
	f2 = range(fmin, step = fmin, fmax)

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
			X₂ = rfft(x₂ .- mean(x₂))
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
		return coherence(f, C̄, P̄₁, P̄₂, n_segments, y₁_err, y₂_err)
	end
	return f, C̄, P̄₁, P̄₂, C
end

@doc raw"""
	coherence(f, C̄, P̄₁, P̄₂, n_segments, y₁_err = nothing, y₂_err = nothing)

Compute the coherence between two time series defined in [1997ApJ...474L..43V](@cite) using the cross-periodogram `C̄`, the mean power spectrum of the first time series `P̄₁`, and the mean power spectrum of the second time series `P̄₂`. The coherence is defined as:

```math
\gamma^2 (f) = \frac{|\overline{C}(f)|^2}{\overline{P}_1(f) \overline{P}_2(f)}
```
with the error in the coherence given by:

```math
\sigma_{\gamma^2} = \sqrt{2} \sqrt{\gamma^2} \frac{1 - \gamma^2}{\sqrt{M}}
```

where `M` is the number of segments over which the coherence is computed.

The phase lag is defined as:

```math
\Delta \phi (f) = \arg \overline{C}(f)
```

with the error in the phase lag given by:

```math
\sigma_{\Delta \phi} = \frac{1 - \gamma^2}{\sqrt{2 \gamma^2 M}}
```

The time lag is defined as:

```math
\Delta \tau (f) = \frac{\Delta \phi}{2 \pi f}
```

with the error in the time lag given by:

```math
\sigma_{\Delta \tau} = \frac{\sigma_{\Delta \phi}}{2 \pi f}
```

If the errors in the time series are provided, the corrected coherence is computed using the method described in [1997ApJ...474L..43V](@cite).

# Arguments
- `f::Array{Float64, 1}`: Frequency array.
- `C̄::Array{Complex{Float64}, 1}`: Cross-periodogram.
- `P̄₁::Array{Float64, 1}`: Mean power spectrum of the first time series.
- `P̄₂::Array{Float64, 1}`: Mean power spectrum of the second time series.
- `n_segments::Int`: Number of segments over which the coherence is computed.
- `y₁_err::Array{Float64, 1}`: Error in the first time series.
- `y₂_err::Array{Float64, 1}`: Error in the second time series.

if the errors in the time series are not provided, the function returns the coherence, phase lag, and time lag. Otherwise, it returns the corrected coherence, phase lag, and time lag.

# Returns
- `f::Array{Float64, 1}`: Frequency array.
- `γ²::Array{Float64, 1}`: Coherence.
- `Δφ::Array{Float64, 1}`: Phase lag.
- `γ²_err::Array{Float64, 1}`: Error in the coherence.
- `Δφ_err::Array{Float64, 1}`: Error in the phase lag.
- `Δτ::Array{Float64, 1}`: Time lag.
- `Δτ_err::Array{Float64, 1}`: Error in the time lag.
- `P̄₁::Array{Float64, 1}`: Mean power spectrum of the first time series.
"""
function coherence(f, C̄, P̄₁, P̄₂, n_segments, y₁_err = nothing, y₂_err = nothing)
	# compute the raw coherence
	γ² = abs.(C̄) .^ 2 ./ (P̄₁ .* P̄₂)
	γ²_err = √2 .* .√abs.(γ²) .* (1.0 .- γ²) ./ √n_segments # Eq. 9.81 in Bendat and Piersol (2010) as square root of the variance

	

	# compute the phase lags
	Δφ = angle.(C̄) # between -π and π
	Δφ_err = (1 .- γ²) ./ .√(2 .* abs.(γ²) .* n_segments)

	# compute the time lags
	Δτ = Δφ ./ (2π .* f)
	Δτ_err = Δφ_err ./ (2π .* f)

	if !isnothing(y₁_err) && !isnothing(y₂_err)
		# compute corrected coherence

		# noise powers
		N₁ = size(y₁_err, 1) * mean(y₁_err .^ 2)
		N₂ = size(y₁_err, 1) * mean(y₂_err .^ 2)

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
		return f, γ², γ²_corrected, Δφ, γ²_err, γ²_corrected_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, N₁, N₂, n²

	end
	return f, γ², Δφ, γ²_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, C̄
end