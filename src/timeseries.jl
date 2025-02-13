using Unitful


abstract type TimeStamps end

"""
	IrregularTimeStamps

A structure to store the time stamps of a time series. It contains the time stamps, the time step, the duration of the time series and a boolean to indicate if the time stamps are irregular.
If the time stamps are irregular, the time step is the minimum time step between two time stamps.

# Fields
- `time::AbstractVector{T}`: the time stamps
- `unit::Unitful.Unit`: the unit of the time stamps
- `timezero::Td`: the time of the first time stamp
- `Δt::T`: the time step
- `duration::T`: the duration of the time series
"""
struct IrregularTimeStamps{T <: Float64, Tu <: Unitful.Units, Tz} <: TimeStamps
    time::AbstractVector{T}
    unit::Tu
    timezero::Tz # can be a date or a time
    Δt::T
    duration::T

    function IrregularTimeStamps(time::AbstractVector{T}, unit::Tu, timezero::Tz, Δt::T, duration::T) where {T <: Float64, Tu <: Unitful.Units, Tz}
        @assert dimension(unit) == Unitful.𝐓 "The unit of the time stamps does not have a Time (Unitful.𝐓) dimension. https://painterqubits.github.io/Unitful.jl/stable/defaultunits/#Base-dimensions"
        return new{T, Tu, Tz}(time, unit, timezero, Δt, duration)
    end

    function IrregularTimeStamps(time::AbstractVector{T}, unit::Tu, timezero::Tz) where {T <: Float64, Tu <: Unitful.Units, Tz}
        @assert dimension(unit) == Unitful.𝐓 "The unit of the time stamps does not have a Time (Unitful.𝐓) dimension. https://painterqubits.github.io/Unitful.jl/stable/defaultunits/#Base-dimensions"
        Δt = minimum(diff(time))
        duration = time[end] - time[1]
        return new{T, Tu, Tz}(time, unit, timezero, Δt, duration)
    end

    function IrregularTimeStamps(time::AbstractVector{T}, unit::Tu) where {T <: Float64, Tu <: Unitful.Units}
        @assert dimension(unit) == Unitful.𝐓 "The unit of the time stamps does not have a Time (Unitful.𝐓) dimension. https://painterqubits.github.io/Unitful.jl/stable/defaultunits/#Base-dimensions"

        Δt = minimum(diff(time))
        duration = time[end] - time[1]
        timezero = time[1]
        return new{T, Tu, Tz}(time, unit, timezero, Δt, duration)
    end
end

"""
	RegularTimeStamps

A structure to store the time stamps of a time series. It contains the time stamps, the time step, the duration of the time series and a boolean to indicate if the time stamps are irregular.
If the time stamps are irregular, the time step is the minimum time step between two time stamps.

# Fields
- `Δt::T`: the time step
- `duration::T`: the duration of the time series
- `unit::Unitful.Unit`: the unit of the time stamps
- `timezero::Td`: the time of the first time stamp
"""
struct RegularTimeStamps{T, Tu <: Unitful.Units, Tz} <: TimeStamps
    Δt::T
    duration::T
    unit::Tu
    timezero::Tz # can be a date or a time

    function RegularTimeStamps(Δt::T, duration::T, unit::Tu, timezero::Tz) where {T, Tu <: Unitful.Units, Tz}
        @assert dimension(unit) == Unitful.𝐓 "The unit of the time stamps does not have a Time (Unitful.𝐓) dimension. See https://painterqubits.github.io/Unitful.jl/stable/defaultunits/#Base-dimensions"
        return new{T, Tu, T}(Δt, duration, unit, timezero)
    end

    function RegularTimeStamps(Δt::T, duration::T, unit::Tu) where {T, Tu <: Unitful.Units}
        @assert dimension(unit) == Unitful.𝐓 "The unit of the time stamps does not have a Time (Unitful.𝐓) dimension. See https://painterqubits.github.io/Unitful.jl/stable/defaultunits/#Base-dimensions"
        timezero = 0.0
        return new{T, Tu, T}(Δt, duration, unit, timezero)
    end
end

"""
	TimeSeriesData

	A type for the data in the time series.

	# Fields
	- `columns::Vector{Symbol}`: the names of the columns
	- `data::AbstractArray{T, 2}`: the data
	- `units::Vector{Unitful.Unit}`: the units of the columns
"""
struct TimeSeriesData{T, Tu <: Unitful.Units}
    columns::Vector{Symbol}
    data::AbstractArray{T, 2}
    units::Vector{Tu}

    function TimeSeriesData(columns::Vector{Symbol}, data::AbstractArray{T, 2}, units::Vector{Tu}) where {T, Tu <: Unitful.Units}
        @assert length(columns) == size(data, 2) "The number of columns in the data does not match the number of columns in the columns vector."
        @assert length(units) == size(data, 2) "The number of columns in the data does not match the number of columns in the units vector."

        # check that we don't have invalid values
        for i in 1:size(data, 2)
            @assert all(.!isnan.(data[:, i])) "The column $columns[i] has NaN values."
            @assert all(.!isinf.(data[:, i])) "The column $columns[i] has Inf values."
        end

        return new{T, Tu}(columns, data, units)
    end
end

"""
	TimeSeries

	A structure to store a time series. It contains the time stamps, the data, and a dictionary to store metadata.

	# Fields
	- `time::TimeStamps{Tt,Tz}`: the time stamps
	- `data::TimeSeriesData{Td, N}`: the data
	- `metadata::Dict{Symbol, Any}`: a dictionary to store metadata

	# Entries in metadata

	- `:name`: a string to store the name of the time series
	- `:description`: a string to store a description of the time series
	- `:unit`: a string to store the unit of the data
	- `:unit_time`: a string to store the unit of the time stamps
	- `:timezero`: a string to store the time of the first time stamp
	- `:fake`: a boolean to indicate if the time series is fake
	- `:instrument`: a string to store the instrument used to measure the time series
	- `:irregular`: a boolean to indicate if the time series has irregular time stamps
"""
struct TimeSeries{T, Typet <: TimeStamps, Tu <: Unitful.Units}
    time::Typet
    data::TimeSeriesData{T, Tu}
    metadata::Dict{Symbol, Any}
end

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
