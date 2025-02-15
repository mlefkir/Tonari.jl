abstract type TimeStamps end

"""
	IrregularTimeStamps

A structure to store the time stamps of a time series. It contains the time stamps, the time step, the duration of the time series and a boolean to indicate if the time stamps are irregular.
If the time stamps are irregular, the time step is the minimum time step between two time stamps.

# Fields
- `time::AbstractVector{T}`: the time stamps
- `unit::Unitful.Unit`: the unit of the time stamps (e.g. seconds, days, etc)
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
        @assert dimension(unit) == Unitful.𝐓 "The unit of the time stamps does not have a Time (Unitful.𝐓) dimension. https://painterqubits.github.io/Unitful.jl/stable/defaultunits/#Base-dimensions" # COV_EXCL_LINE
        Δt = minimum(diff(time))
        duration = time[end] - time[1]
        return new{T, Tu, Tz}(time, unit, timezero, Δt, duration)
    end

    function IrregularTimeStamps(time::AbstractVector{T}, unit::Tu) where {T <: Float64, Tu <: Unitful.Units}
        @assert dimension(unit) == Unitful.𝐓 "The unit of the time stamps does not have a Time (Unitful.𝐓) dimension. https://painterqubits.github.io/Unitful.jl/stable/defaultunits/#Base-dimensions" # COV_EXCL_LINE
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
- `unit::Unitful.Unit`: the unit of the time stamps (e.g. seconds, days, etc)
- `timezero::Td`: the time of the first time stamp
"""
struct RegularTimeStamps{T, Tu <: Unitful.Units, Tz} <: TimeStamps
    Δt::T
    duration::T
    unit::Tu
    timezero::Tz # can be a date or a time

    function RegularTimeStamps(Δt::T, duration::T, unit::Tu, timezero::Tz) where {T, Tu <: Unitful.Units, Tz}
        @assert dimension(unit) == Unitful.𝐓 "The unit of the time stamps does not have a Time (Unitful.𝐓) dimension. See https://painterqubits.github.io/Unitful.jl/stable/defaultunits/#Base-dimensions" # COV_EXCL_LINE
        return new{T, Tu, T}(Δt, duration, unit, timezero)
    end

    function RegularTimeStamps(Δt::T, duration::T, unit::Tu) where {T, Tu <: Unitful.Units}
        @assert dimension(unit) == Unitful.𝐓 "The unit of the time stamps does not have a Time (Unitful.𝐓) dimension. See https://painterqubits.github.io/Unitful.jl/stable/defaultunits/#Base-dimensions" # COV_EXCL_LINE
        timezero = 0.0
        return new{T, Tu, T}(Δt, duration, unit, timezero)
    end
end

"""
	TimeSeriesData

	A type for the data in the time series.

	# Fields
	- `Union{Vector{Symbol},Symbol}`: the names of the columns, e.g. flux, error. Can be a vector of symbols or a single symbol if there is only one column.
	- `data::Union{AbstractArray{T,2},Vector{T}}`: the data array. Can be a matrix of data or a vector if there is only one column.
	- `units::Union{Vector{Tu},Tu}`: the units of the data. Can be a vector of units or a single unit if there is only one column or the same unit is used for all columns.
"""
struct TimeSeriesData{T, Tu <: Unitful.Units}
    columns::Union{Vector{Symbol}, Symbol}
    data::Union{AbstractArray{T, 2}, Vector{T}}
    units::Union{Vector{Tu}, Tu}

    function TimeSeriesData(columns::Union{Vector{Symbol}, Symbol}, data::Union{AbstractArray{T, 2}, Vector{T}}, units::Union{Vector{Tu}, Tu}) where {T, Tu <: Unitful.Units}

        if columns isa Vector
            @assert length(columns) == size(data, 2) "The number of columns in the data does not match the number of columns in the columns vector."
        end
        if units isa Vector
            @assert length(units) == size(data, 2) "The number of columns in the data does not match the number of columns in the units vector."
        end
        if data isa Vector{T}
            @assert all(.!isnan.(data)) "The column $columns[i] has NaN values."
            @assert all(.!isinf.(data)) "The column $columns[i] has Inf values."
        else
            @assert columns isa Vector && length(columns) == size(data, 2) "The number of columns in the data does not match the number of columns in the columns vector."
            # check that we don't have invalid values
            for i in eachindex(view(data, 1, :))
                @assert all(.!isnan.(data[:, i])) "The column $columns[i] has NaN values."
                @assert all(.!isinf.(data[:, i])) "The column $columns[i] has Inf values."
            end
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
	- `:unit`: a string to store the unit(s) of the data
	- `:unit_time`: unit of the time stamps (default: seconds)
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

"""
    TimeSeries(time::AbstractVector{T}, data::AbstractVector{Td}, metadata::Dict{Symbol,Any} = Dict{Symbol,Any}(), atol::T=0) where {T, Td}

    A constructor for the TimeSeries structure.

    This function lets you create a TimeSeries structure using the times, data, and optional metadata.

    # Arguments
    - `time::AbstractVector{T}`: the time stamps of the time series (length N)
    - `data::AbstractVector{Td}`: the data array (N x M) where N is the number of time stamps and M is the number of columns, M can be 1 for a single column, then data is a vector
    - `metadata::Dict{Symbol,Any}`: a dictionary to store metadata
    - `atol::T=0`: the absolute tolerance for the time stamps to check if they are regular

    # Metadata (optional but highly recommended)
    - `:irregular`::Bool : a boolean to indicate if the time stamps are irregular
    - `:unit_time`: the unit of the time stamps (default: seconds) using Unitful
    - `:timezero`: the time of the first time stamp (default: time[1])
    - `:columns`: the names of the columns (length M, default: y1, y2, ...), if M = 1, then the column is a Symbol
    - `:unit`: the units of the columns (using Unitful) (length M, default: dimensionless), if M = 1, then the unit is a Symbol, but if all columns have the same unit, then the unit is a single unit
    - `:name`: the name of the time series (optional)
    - `:description`: a description of the time series (optional)
    - `:fake`: a boolean to indicate if the time series is fake (optional)
    - `:instrument`: the instrument used to measure the time series (optional)

    # Returns
    - `TimeSeries`: a TimeSeries structure with the time stamps, data, and metadata
"""
function TimeSeries(
        time::AbstractVector{T},
        data::Union{AbstractArray{Td, 2}, Vector{Td}},
        metadata::Dict{Symbol, Any} = Dict{Symbol, Any}();
        atol::T = 0.0
    ) where {T, Td}

    ### sanity checks
    @assert length(time) == size(data, 1) "The time stamps and the data do not have the same length."
    @assert issorted(time) "The time stamps are not sorted in ascending order."

    if metadata == Dict{Symbol, Any}()
        @info "No metadata provided. Creating an empty metadata dictionary."
        metadata = Dict{Symbol, Any}()
    end
    ### work on the time stamps
    # check if the time stamps are regular
    regular_sampling = regular_time_sampling(time, atol = atol)
    if :irregular in keys(metadata)
        @assert metadata[:irregular] == !regular_sampling "The metadata and the time stamps disagree on the time sampling.
        metadata:$(metadata[:irregular]) != regular_sampling $(!regular_sampling)"
    else
        metadata[:irregular] = !regular_sampling
    end

    # check if units for the time stamps are provided
    if :unit_time in keys(metadata)
        unit_time = metadata[:unit_time]
    else
        @info "No unit for the time stamps `unit_time` provided. Assuming seconds."
        unit_time = u"s"
        metadata[:unit_time] = unit_time
    end
    # check if the timezero is provided
    if :timezero in keys(metadata)
        timezero = metadata[:timezero]
    else
        @info "No timezero provided. Assuming the first time stamp is the timezero."
        timezero = time[1]
        metadata[:timezero] = timezero
    end
    # build the time stamps
    if regular_sampling
        Δt = mean(diff(time))
        duration = time[end] - time[1]
        time_stamps = RegularTimeStamps(Δt, duration, unit_time, timezero)
    else
        time_stamps = IrregularTimeStamps(time, unit_time, timezero)
    end

    ### work on the data
    if :columns in keys(metadata)
        columns = metadata[:columns]
    else
        @info "No columns provided. Assuming the columns are y1, y2, ..."
        columns = [Symbol("y$i") for i in eachindex(view(data, 1, :))]
        metadata[:columns] = columns
    end

    if :units in keys(metadata)
        units = metadata[:units]
    else
        @info "No units for the data provided. Assuming dimensionless."
        units = [NoUnits for i in eachindex(view(data, 1, :))]
        metadata[:units] = units
    end

    # make the TimeSeriesData
    ts_data = TimeSeriesData(columns, data, units)
    return TimeSeries(time_stamps, ts_data, metadata)
end

@doc raw"""
	fill_gaps(rng::AbstractRNG, t::Vector{Float64}, x::Vector{Float64}, σₓ = nothing; randomise_values::Bool = true, poisson::Bool = true)

    Fill gaps in a time series data with linear interpolation. If `randomise_values = true`, the interpolated values are drawn from a normal distribution with the mean and standard deviation of the data. If `poisson = true`, the interpolated values are drawn from a Poisson distribution with the mean of the data.

    # Arguments
    - `rng::AbstractRNG`: random number generator
    - `t::Vector{Float64}`: time array
    - `x::Vector{Float64}`: data array
    - `σₓ`::Vector{Float64}: standard deviation of the data array (optional)
    - `randomise_values::Bool`: whether to randomise the interpolated values (default: true)
    - `poisson::Bool`: whether to use a Poisson distribution for the interpolated values, this assumes that the data is in counts/dt units (default: true)
    - `Δt::Float64`: the time step of the time series (optional)

    # Returns
    - `t_filled::Vector{Float64}`: the time array with the gaps filled
    - `x_filled::Vector{Float64}`: the data array with the gaps filled
    - `σ_filled::Vector{Float64}`: the standard deviation of the data array with the gaps filled (optional)
"""
function fill_gaps(
        rng::AbstractRNG,
        t::AbstractVector{Float64},
        x::AbstractVector{Float64},
        σₓ = nothing;
        randomise_values::Bool = true,
        poisson::Bool = true,
        Δt = nothing
    )

    # check the time series
    time_series_sanity_checks(t, x, σₓ)

    Δt = isnothing(Δt) ? median(diff(t)) : Δt
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
    if !isnothing(σₓ) && randomise_values
        σ_filled = vcat(σₓ, σ_interp)
        σ_filled = σ_filled[p]
        return t_filled, x_filled, σ_filled
    end

    return t_filled, x_filled
end

fill_gaps(t::AbstractVector{Float64}, x::AbstractVector{Float64}, σₓ = nothing; randomise_values::Bool = true, poisson::Bool = true) = fill_gaps(Random.GLOBAL_RNG, t, x, σₓ; randomise_values = randomise_values, poisson = poisson)

"""
	time_series_sanity_checks(t₁::AbstractVector{T}, x₁::AbstractVector{T}, σ₁=nothing) where T

	Perform sanity checks on the time series. The function checks the following:

	- The time and value vectors have the same length
	- The time and uncertainty vectors have the same length
	- The time series are sorted in ascending order
	- The time series do not contain infinities or NaNs

	# Arguments
	- `t₁::AbstractVector{T}`: time points of the time series
	- `x₁::AbstractVector{T}`: values of the time series
	- `σ₁::AbstractVector{T}=nothing`: uncertainty of the time series
"""
function time_series_sanity_checks(
        t₁::AbstractVector{T},
        x₁::AbstractVector{T},
        σ₁ = nothing
    ) where {T}

    # assert that the time series have the same length
    @assert length(t₁) == length(x₁) "Time and value vectors of the time series must have the same length"
    if !isnothing(σ₁)
        @assert length(t₁) == length(σ₁) "Time and uncertainty vectors of the time series must have the same length"
    end

    # check the order of the time series
    if !issorted(t₁)
        @info "Sorting the time series"
        p = sortperm(t₁)
        t₁ = t₁[p]
        x₁ = x₁[p]
        if !isnothing(σ₁)
            σ₁ = σ₁[p]
        end
    end
    # check for infinities and NaNs
    @assert !(any(isinf.(x₁)) || any(isnan.(x₁))) "The time series contains infinities or NaNs"

    return if !isnothing(σ₁)
        @assert !(any(isinf.(σ₁)) || any(isnan.(σ₁))) "Uncertainty on the time series contains infinities or NaNs"
    end
end
