using StatsBase

"""
	time_series_sanity_checks(t₁::AbstractVector{T}, x₁::AbstractVector{T}, σ₁=nothing) where T

	Perform sanity checks on the time series. The function checks the following:

	- The time and value vectors have the same length
	- The time and uncertainty vectors have the same length
	- The time series are sorted in ascending order
	- The time series do not contain infinities or NaNs

	# Arguments
	- `t₁::AbstractVector{T}`: time points of the first time series
	- `x₁::AbstractVector{T}`: values of the first time series
	- `σ₁::AbstractVector{T}=nothing`: uncertainty of the first time series
"""
function time_series_sanity_checks(
        t₁::AbstractVector{T},
        x₁::AbstractVector{T},
        σ₁ = nothing
    ) where {T}

    # assert that the time series have the same length
    @assert length(t₁) == length(x₁) "Time and value vectors of the first time series must have the same length"
    if !isnothing(σ₁)
        @assert length(t₁) == length(σ₁) "Time and uncertainty vectors of the first time series must have the same length"
    end

    # check the order of the time series
    if !issorted(t₁)
        @info "Sorting the first time series"
        p = sortperm(t₁)
        t₁ = t₁[p]
        x₁ = x₁[p]
        if !isnothing(σ₁)
            σ₁ = σ₁[p]
        end
    end
    # check for infinities and NaNs
    if any(isinf.(x₁)) || any(isnan.(x₁))
        @error "First time series contains infinities or NaNs"
    end

    return if !isnothing(σ₁)
        if any(isinf.(σ₁)) || any(isnan.(σ₁))
            @error "Uncertainty of the first time series contains infinities or NaNs"
        end
    end
end

"""
	cross_correlate(t₁::AbstractVector{T}, x₁::AbstractVector{T}, t₂::AbstractVector{T}, x₂::AbstractVector{T}; σ₁=nothing, σ₂=nothing, τ_list=nothing, Δτ=nothing, max_lag=nothing, local_estimate=false, both_ways=true, method="iccf") where T

	Compute the cross-correlation function between two time series.

	The function computes the cross-correlation function between two time series by linearly interpolating the second time series using the interpolated cross-correlation function (ICCF) method.
	Most of this code is based on the R package `sour` by Simon Vaughan (University of Leicester).



	# Arguments
	- `t₁::AbstractVector{T}`: time points of the first time series
	- `x₁::AbstractVector{T}`: values of the first time series
	- `t₂::AbstractVector{T}`: time points of the second time series
	- `x₂::AbstractVector{T}`: values of the second time series

	# Optional arguments
	- `σ₁::AbstractVector{T} = nothing`: uncertainty of the first time series
	- `σ₂::AbstractVector{T} = nothing`: uncertainty of the second time series
	- `τ_list::AbstractVector{T} = nothing`: list of time lags
	- `Δτ::T = nothing`: time step between the time lags
    - `max_lag::T = nothing`: maximum time lag
    - `local_estimate::Bool = false`: flag to estimate the local mean and standard deviation of the time series
    - `both_ways::Bool = true`: flag to compute the cross-correlation function in both ways
    - `method::String = "iccf"`: method to compute the cross-correlation function
"""
function cross_correlate(
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
        method = "iccf"
    ) where {T}

    # sanity checks
    time_series_sanity_checks(t₁, x₁, σ₁)
    time_series_sanity_checks(t₂, x₂, σ₂)

    # if the time lags are not provided, compute them
    if !isnothing(τ_list)
        n_lags = length(τ_list)
        # Δτ = τ_list[2]-τ_list[1]
        # max_lag = (length(τ_list)-1)*Δτ/2
    else
        Δτ = isnothing(Δτ) ? minimum(diff(t₁)) : Δτ
        max_lag = isnothing(max_lag) ? (maximum([t₁; t₂]) - minimum([t₁; t₂])) / 4 : max_lag

        if max_lag > (maximum([t₁; t₂]) - minimum([t₁; t₂])) / 4
            @warn "Max lag is larger than the time range of the time series setting it to the maximum time range divided by 4"
            max_lag = maximum([t₁; t₂]) - minimum([t₁; t₂]) / 4
        end

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


    return τ_list, r
end


"""
	iccf(t₁::AbstractVector{T}, x₁::AbstractVector{T}, t₂::AbstractVector{T}, x₂::AbstractVector{T}, τ_list::AbstractVector{T} ; local_estimate=false) where T

	Compute the interpolated cross-correlation function (ICCF) between two time series.


	The function computes the cross-correlation function between two time series by linearly interpolating the second time series
	on the time points of the first time series. The ICCF is presented in Edelson et al. 2017  https://ui.adsabs.harvard.edu/abs/2017ApJ...840...41E/abstract
    and Peterson et al. 2004 https://ui.adsabs.harvard.edu/abs/2004ApJ...613..682P/abstract.
    The current implementation is based on the R package `sour` by Simon Vaughan (University of Leicester).

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
