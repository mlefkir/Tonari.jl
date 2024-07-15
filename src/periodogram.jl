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
