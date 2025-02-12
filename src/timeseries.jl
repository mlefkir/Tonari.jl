@doc raw"""
	fill_gaps(rng::AbstractRNG, t::Vector{Float64}, x::Vector{Float64}, σₓ = nothing; randomise_values::Bool = true, poisson::Bool = true)

Fill gaps in a time series data with linear interpolation. If `randomise_values = true`, the interpolated values are drawn from a normal distribution with the mean and standard deviation of the data. If `poisson = true`, the interpolated values are drawn from a Poisson distribution with the mean of the data.

# Arguments
- `rng::AbstractRNG`: random number generator
- `t::Vector{Float64}`: time array
- `x::Vector{Float64}`: data array
- `σₓ`::Vector{Float64}: standard deviation of the data array (optional)
- `randomise_values::Bool`: whether to randomise the interpolated values
- `poisson::Bool`: whether to use a Poisson distribution for the interpolated values, this assumes that the data is in counts/dt
"""
function fill_gaps(rng::AbstractRNG, t::AbstractVector{Float64}, x::AbstractVector{Float64}, σₓ = nothing; randomise_values::Bool = true, poisson::Bool = true)

    Δt = minimum(unique(diff(t))) # minimum time step
    @info "Filling gaps assuming Δt = $Δt"
    # time array we should have
    t_new = range(t[1], stop = t[end], step = Δt)
    # find times we are missing in t[i]
    missing_times = setdiff(t_new, t)
    @assert size(missing_times, 1) > 0 "No missing times found check the time array"
    # get the interpolated values
    x_interp = linear_interpolation(t, x)(missing_times)
    # randomise the interpolated values
    if randomise_values
        if poisson
            x_interp = rand.(rng, Poisson.(x_interp * Δt)) / Δt
            if !isnothing(σₓ)
                σ_interp = sqrt.(x_interp * Δt) / Δt
            end
        else
            error("Not implemented")
        end
    end


    # append to the original data
    x_filled = vcat(x, x_interp)
    t_filled = vcat(t, missing_times)

    # sort
    p = sortperm(t_filled)
    t_filled = t_filled[p]
    x_filled = x_filled[p]
    if !isnothing(σₓ)
        σ_filled = vcat(σₓ, σ_interp)
        σ_filled = σ_filled[p]
        return t_filled, x_filled, σ_filled
    end

    return t_filled, x_filled
end

fill_gaps(t::AbstractVector{Float64}, x::AbstractVector{Float64}, σₓ = nothing; randomise_values::Bool = true, poisson::Bool = true) = fill_gaps(Random.GLOBAL_RNG, t, x, σₓ; randomise_values = randomise_values, poisson = poisson)
