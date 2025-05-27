"""
    get_lag_list(Δτ, max_lag,t₁,t₂,τ_list)

Get the list of time lags to compute the cross-correlation function.

# Arguments
- `Δτ::T`: time step between the time lags
- `max_lag::T`: maximum time lag
- `t₁::AbstractVector{T}`: time stamps of the first time series
- `t₂::AbstractVector{T}`: time stamps of the second time series
- `τ_list::AbstractVector{T}`: list of time lags

# Returns
- `τ_list::AbstractVector{T}`: list of time lags
- `Δτ::T`: time step between the time lags
- `max_lag::T`: maximum time lag
- `n_lags::Int64`: number of time lags
"""
function get_lag_list(Δτ, max_lag, t₁, t₂, τ_list)
    # if the time lags are not provided, compute them
    if !isnothing(τ_list)
        n_lags = length(τ_list)
        # Δτ = τ_list[2]-τ_list[1]
        # max_lag = (length(τ_list)-1)*Δτ/2
    else
        Δτ = isnothing(Δτ) ? minimum(diff(t₁)) : Δτ
        max_lag = isnothing(max_lag) ? (maximum([t₁; t₂]) - minimum([t₁; t₂])) / 4 : max_lag

        # if max_lag > (maximum([t₁; t₂]) - minimum([t₁; t₂])) / 4
        #     @warn "Max lag is larger than the time range of the time series setting it to the maximum time range divided by 4"
        #     max_lag = maximum([t₁; t₂]) - minimum([t₁; t₂]) / 4
        # end

        # compute the number of lags
        n_lags = ceil(Int64, 2 * max_lag / Δτ)
        # make sure n_lags is odd, so that tau=0 is included
        if n_lags % 2 == 0
            n_lags += 1
        end
        # update the max_lag
        max_lag = (n_lags - 1) * Δτ / 2

        # create the tau array
        τ_list = -max_lag:Δτ:max_lag
    end
    return τ_list, Δτ, max_lag, n_lags
end

"""
	cross_correlate(t₁::AbstractVector{T}, x₁::AbstractVector{T}, t₂::AbstractVector{T}, x₂::AbstractVector{T}; σ₁=nothing, σ₂=nothing, τ_list=nothing, Δτ=nothing, max_lag=nothing, local_estimate=false, both_ways=true, method="iccf", compute_errors=false, peak_frac=0.8, bootstrap=false, n_simulations=1_000,
    skip_sanity_checks = false) where T

Compute the cross-correlation function between two time series.

The function computes the cross-correlation function between two time series by linearly interpolating the second time series using the interpolated cross-correlation function (ICCF) method.
Most of this code is based on the R package `sour` by Simon Vaughan, University of Leicester https://github.com/SimonVaughanDataAndCode/sour.

# Arguments
- `t₁::AbstractVector{T}`: time stamps of the first time series
- `x₁::AbstractVector{T}`: values of the first time series
- `t₂::AbstractVector{T}`: time stamps of the second time series
- `x₂::AbstractVector{T}`: values of the second time series

# Optional arguments
- `σ₁::AbstractVector{T} = nothing`: uncertainty of the first time series
- `σ₂::AbstractVector{T} = nothing`: uncertainty of the second time series
- `τ_list::AbstractVector{T} = nothing`: list of time lags, if not provided, it will be computed based on the time range of the time series
- `Δτ::T = nothing`: time step between the time lags
- `max_lag::T = nothing`: maximum time lag
- `local_estimate::Bool = false`: flag to estimate the local mean and standard deviation of the time series
- `both_ways::Bool = true`: Compute the cross-correlation function in both ways, i.e., from the first to the second time series and from the second to the first time series
- `method::String = "iccf"`: method to compute the cross-correlation function
- `compute_errors::Bool = false`: flag to compute the errors of the cross-correlation function
- `peak_frac::T = 0.8`: fraction of the peak to integrate to get the lag
- `bootstrap::Bool = false`: flag to compute the errors using bootstrapping
- `n_simulations::Int = 1_000`: number of simulations to compute the errors
- `skip_sanity_checks::Bool = false`: flag to skip the sanity checks

# Returns
- `τ_list::AbstractVector{T}`: list of time lags
- `r::AbstractVector{T}`: cross-correlation function
- `q::AbstractVector{T}`: Distribution of centroid lags if `compute_errors` is `true`
"""
function cross_correlate(
        rng::AbstractRNG,
        t₁::AbstractVector{T},
        x₁::AbstractVector{T},
        t₂::AbstractVector{T},
        x₂::AbstractVector{T};
        σ₁ = nothing,
        σ₂ = nothing,
        τ_list = nothing,
        Δτ = nothing,
        max_lag = nothing,
        local_estimate = false,
        both_ways = true,
        method = "iccf",
        compute_errors = false,
        peak_frac = 0.8,
        bootstrap = false,
        n_simulations = 1_000,
        skip_sanity_checks = false
    ) where {T}

    # sanity checks
    if !skip_sanity_checks
        time_series_sanity_checks(t₁, x₁, σ₁)
        time_series_sanity_checks(t₂, x₂, σ₂)
    end

    # get array of lag times
    τ_list, Δτ, max_lag, n_lags = get_lag_list(Δτ, max_lag, t₁, t₂, τ_list)

    # interpolate the second time series
    if method == "iccf"
        if both_ways
            r₁, _ = iccf(t₁, x₁, t₂, x₂, τ_list, local_estimate = local_estimate)
            r₂, _ = iccf(t₂, x₂, t₁, x₁, -τ_list, local_estimate = local_estimate)
            r = (r₁ + r₂) / 2
        else
            r, _ = iccf(t₁, x₁, t₂, x₂, τ_list, local_estimate = local_estimate)
        end
    else
        error("Method not implemented")
    end

    # sanity check
    if any(r .> 1) || any(r .< -1)
        @warn "Cross-correlation function is not within [-1,1]"
    end

    # compute the errors
    if compute_errors
        q = iccf_errors(
            rng, t₁, x₁, t₂, x₂, σ₁ = σ₁, σ₂ = σ₂,
            peak_frac = peak_frac,
            local_estimate = local_estimate,
            τ_list = τ_list,
            bootstrap = bootstrap,
            n_simulations = n_simulations
        )

        return τ_list, r, q

    end
    return τ_list, r
end


cross_correlate(
    t₁::AbstractVector,
    x₁::AbstractVector,
    t₂::AbstractVector,
    x₂::AbstractVector;
    σ₁ = nothing,
    σ₂ = nothing,
    τ_list = nothing,
    Δτ = nothing,
    max_lag = nothing,
    local_estimate = false,
    both_ways = true,
    method = "iccf",
    compute_errors = false,
    peak_frac = 0.8,
    bootstrap = false,
    n_simulations = 1_000,
    skip_sanity_checks = false
) = cross_correlate(
    Random.GLOBAL_RNG,
    t₁,
    x₁,
    t₂,
    x₂;
    σ₁ = σ₁,
    σ₂ = σ₂,
    τ_list = τ_list,
    Δτ = Δτ,
    max_lag = max_lag,
    local_estimate = local_estimate,
    both_ways = both_ways,
    method = method,
    compute_errors = compute_errors,
    peak_frac = peak_frac,
    bootstrap = bootstrap,
    n_simulations = n_simulations,
    skip_sanity_checks = skip_sanity_checks
)


"""
	iccf(t₁::AbstractVector{T}, x₁::AbstractVector{T}, t₂::AbstractVector{T}, x₂::AbstractVector{T}, τ_list::AbstractVector{T} ; local_estimate=false) where T

Compute the interpolated cross-correlation function (ICCF) between two time series.

The function computes the cross-correlation function between two time series by linearly interpolating the second time series
on the time points of the first time series. The ICCF is presented in Edelson et al. 2017  https://ui.adsabs.harvard.edu/abs/2017ApJ...840...41E/abstract
and Peterson et al. 2004 https://ui.adsabs.harvard.edu/abs/2004ApJ...613..682P/abstract.
The current implementation is based on the R package `sour` by Simon Vaughan, University of Leicester https://github.com/SimonVaughanDataAndCode/sour.

https://ui.adsabs.harvard.edu/abs/1986ApJ...305..175G/abstract

# Arguments
- `t₁::AbstractVector{T}`: time points of the first time series
- `x₁::AbstractVector{T}`: values of the first time series
- `t₂::AbstractVector{T}`: time points of the second time series
- `x₂::AbstractVector{T}`: values of the second time series
- `τ_list::AbstractVector{T}`: list of time lags
- `local_estimate::Bool=false`: flag to estimate the local mean and standard deviation of the time series

# Returns
- `CCF::Vector{T}`: cross-correlation function
- `nx::Vector{Int64}`: number of points in each time lag bin

# Example
```julia
using Random
Random.seed!(123)
t₁ = 0:0.1:10
x₁ = randn(length(t₁))
t₂ = 0:0.1:10
x₂ = randn(length(t₂))
τ_list = -10:0.1:10
CCF,nx = iccf(t₁, x₁, t₂, x₂, τ_list)
```
"""
function iccf(
        t₁::AbstractVector{T},
        x₁::AbstractVector{T},
        t₂::AbstractVector{T},
        x₂::AbstractVector{T},
        τ_list::AbstractVector{T};
        local_estimate = false
    ) where {T}

    # number of time lags
    n = length(τ_list)
    # preallocate memory for the cross-correlation function
    CCF = Vector{T}(undef, n)
    nx = Vector{Int64}(undef, n)

    # index of the zero lag
    zeroth_lag = trunc(Int64, (length(τ_list) - 1) / 2) + 1

    # linear interpolation of the second time series
    interp_linear = linear_interpolation(t₂, x₂)

    # global mean and standard deviation of the time series
    μ₁, μ₂ = mean(x₁), mean(x₂)
    σ₁, σ₂ = std(x₁), std(x₂)

    for (i, τ) in enumerate(τ_list)

        # lag the first time series
        l = t₁ .- τ

        # mask to select the points within the time range of the second time series
        mask = ((l .> t₂[1]) .&& (l .< t₂[end]))

        # filter the time series
        x₁_m = x₁[mask]

        # interpolate the second time series on the time points of the first time series
        x₂_interp = interp_linear(l[mask])

        # local mean and standard deviation of the time series
        if local_estimate
            μ₁, μ₂ = mean(x₁_m), mean(x₂_interp)
            σ₁, σ₂ = std(x₁_m), std(x₂_interp)
        end

        # number of points in filtered time series
        nᵢ = length(x₁_m)

        # cross-correlation coefficient
        rᵢ = sum((x₁_m .- μ₁) .* (x₂_interp .- μ₂)) ./ (σ₁ * σ₂)

        # store the cross-correlation coefficient and the number of points
        CCF[i] = rᵢ
        nx[i] = nᵢ
    end

    # normalize the cross-correlation function
    if local_estimate
        # by the number of points in each lag bin
        CCF ./= (nx .- 1)
    else
        # by the number of points in the zeroth lag bin
        CCF ./= (nx[zeroth_lag] - 1)
    end

    return CCF, nx
end

"""
    iccf_errors(rng::AbstractRNG, t₁::AbstractVector{T}, x₁::AbstractVector{T}, t₂::AbstractVector{T}, x₂::AbstractVector{T}, σ₁::AbstractVector{T} = nothing, σ₂::AbstractVector{T} = nothing,peak_frac = 0.8,local_estimate = false,τ_list = nothing,bootstrap = false,n_simulations = 1_000) where {T}

Compute the errors of the interpolated cross-correlation function (ICCF) between two time series.

# Arguments
- `rng::AbstractRNG`: Random number generator
- `t₁::AbstractVector{T}`: time points of the first time series
- `x₁::AbstractVector{T}`: values of the first time series
- `t₂::AbstractVector{T}`: time points of the second time series
- `x₂::AbstractVector{T}`: values of the second time series
- `σ₁::AbstractVector{T} = nothing`: uncertainty of the first time series
- `σ₂::AbstractVector{T} = nothing`: uncertainty of the second time series
- `peak_frac::T = 0.8`: fraction of the peak to integrate to get the lag
- `local_estimate::Bool = false`: flag to estimate the local mean and standard deviation of the time series
- `τ_list::AbstractVector{T} = nothing`: list of time lags
- `bootstrap::Bool = false`: flag to compute the errors using bootstrapping
- `n_simulations::Int = 1_000`: number of simulations to compute the errors

# Returns
- `q::Vector{T}`: array of lags
"""
function iccf_errors(
        rng::AbstractRNG,
        t₁::AbstractVector{T},
        x₁::AbstractVector{T},
        t₂::AbstractVector{T},
        x₂::AbstractVector{T};
        σ₁::AbstractVector{T} = nothing,
        σ₂::AbstractVector{T} = nothing,
        peak_frac = 0.8,
        local_estimate = false,
        τ_list = nothing,
        bootstrap = false,
        n_simulations = 1_000
    ) where {T}

    # vector of time lags
    q = Vector{T}(undef, n_simulations)

    @showprogress for i in 1:n_simulations
        # resample the time series and randomize the flux
        t_samp, y_samp, _ = randomise_lc_flux(rng, t₁, x₁, σ = σ₁, bootstrap = bootstrap)
        t_samp_2, y_samp_2, _ = randomise_lc_flux(rng, t₂, x₂, σ = σ₂, bootstrap = bootstrap)

        # get array of lag times
        # τ_list, Δτ, max_lag, n_lags = get_lag_list(Δτ, max_lag, t₁, t₂, τ_list)
        # compute the cross-correlation function
        r₊, _ = iccf(t_samp, y_samp, t_samp_2, y_samp_2, τ_list, local_estimate = local_estimate)
        r₋, _ = iccf(t_samp_2, y_samp_2, t_samp, y_samp, -τ_list, local_estimate = local_estimate)
        r = (r₊ + r₋) / 2
        # find where the peak is
        m = r .>= peak_frac * maximum(r)
        # integrate the peak to get the lag
        τ_peak = sum(r[m] .* τ_list[m]) / sum(r[m])
        # store the lag
        q[i] = τ_peak
    end

    return q
end

"""
	randomise_lc_flux(rng::AbstractRNG, t::Vector{Float64}, x::Vector{Float64}; σ::Vector{Float64} = nothing, bootstrap::Bool = false)

Randomise the light curve by sampling with replacement. If `bootstrap` is `true`, the errors are not added to the light curve.

# Arguments
- `rng::AbstractRNG`: Random number generator
- `t::Vector{Float64}`: Time array
- `x::Vector{Float64}`: Flux array
- `σ::Vector{Float64}`: Error array
- `bootstrap::Bool`: If `true`, the errors are not added to the light curve

# Returns
- `t_samp::Vector{Float64}`: Sampled time array
- `x_samp::Vector{Float64}`: Sampled flux array
- `σ_samp::Vector{Float64}`: Sampled error array
"""
function randomise_lc_flux(rng::AbstractRNG, t::AbstractVector{T}, x::AbstractVector{T}; σ::AbstractVector{T} = nothing, bootstrap::Bool = false) where {T}

    # number of points in time series - and also the number of points to sample
    n_points = length(t)
    # indices of points to sample
    idx = sample(rng, 1:n_points, n_points, replace = true)
    # collect the counts of each index
    counts = collect(values(countmap(idx)))
    idx = collect(keys(countmap(idx)))
    # sort the indices
    p = sortperm(idx)
    idx = idx[p]
    counts = counts[p]

    # sample the time series with replacement
    t_samp, x_samp = t[idx], x[idx]

    σ_samp = 0.0
    # if there are errors, sample them too and rescale by sqrt(n)
    if isnothing(σ)
        σ_samp = σ[idx] / sqrt(counts)
    end

    if !bootstrap && !isnothing(σ)
        x_samp .+= σ_samp .* randn(rng, length(x_samp))
    end


    if isnothing(σ)
        return t_samp, x_samp
    else
        return t_samp, x_samp, σ_samp
    end
end
