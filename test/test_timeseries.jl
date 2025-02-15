using Test, Tonari, Random, Unitful

"""
    setup_timeseries()

    Set up the time series for testing generates
    random time series data for testing
"""
function setup_timeseries()
    rng = MersenneTwister(123)
    N₁, N₂ = 502, 140
    T₁, T₂ = 404.4, 332.2
    Δt₁, Δt₂ = 0.52, 0.38

    # sample data
    t₁ = rand(rng, 0:Δt₁:T₁, N₁)
    t₂ = rand(rng, 0:Δt₂:T₂, N₂)
    t₁_sort = sort(t₁)

    # create time series
    x₁ = randn(rng, N₁)
    x₂ = randn(rng, N₂)
    σ₁ = rand(rng, N₁)
    σ₂ = rand(rng, N₂)
    return t₁, x₁, σ₁, t₂, x₂, σ₂, t₁_sort
end

function generate_regulartimeseries()
    rng = MersenneTwister(123)
    T = 404.4
    t = 0:0.52:T
    n = length(t)
    x = randn(rng, n)
    σ = rand(rng, n)
    return t, x, σ
end

"""
    test_sanity_check()

    Test the time_series_sanity_checks function
"""
function test_sanity_check()
    t₁, x₁, σ₁, t₂, x₂, σ₂, t₁_sort = setup_timeseries()

    # add infinities and NaNs
    x₁[2] = Inf
    x₁[5] = NaN
    σ₂[5] = NaN

    # tests
    @test_throws "Time and value vectors of the time series must have the same length" time_series_sanity_checks(t₁_sort, x₂)
    @test_throws "Time and uncertainty vectors of the time series must have the same length" time_series_sanity_checks(t₂, x₂, σ₁)

    @test_logs (:info, "Sorting the time series") time_series_sanity_checks(t₂, x₂)
    @test_throws "The time series contains infinities or NaNs" time_series_sanity_checks(t₁_sort, x₁)
    return @test_throws "Uncertainty on the time series contains infinities or NaNs" time_series_sanity_checks(sort(t₂), x₂, σ₂)
end


function test_RegularTimeSeries()
    t, x, σ = generate_regulartimeseries()

    # initialising the RegularTimeSeries with just data and times
    ts = TimeSeries(t, [x σ])
    # check the types of the data
    @test ts.time isa RegularTimeStamps
    @test ts.data isa TimeSeriesData
    @test ts.metadata isa Dict

    # check the values in the RegularTimeSeries
    @test ts.time.Δt == 0.52
    @test ts.time.duration == t[end]
    @test ts.time.unit == u"s"
    @test ts.time.timezero == 0.0

    # check the values in the TimeSeriesData
    @test ts.data.data == [x σ]
    @test ts.data.columns == [:y1, :y2]
    @test ts.data.units == [NoUnits, NoUnits]

    # check the values in the metadata
    @test ts.metadata[:irregular] == false
    @test ts.metadata[:unit_time] == u"s"
    @test ts.metadata[:timezero] == 0.0
    @test ts.metadata[:columns] == [:y1, :y2]
    return @test ts.metadata[:units] == [NoUnits, NoUnits]
end

function test_RegularTimeSeries_metadata()
    t, x, σ = generate_regulartimeseries()

    metadata = Dict(
        :name => "Time series",
        :description => "A time series for testing",
        :units => u"erg/cm^2/s",
        :columns => [:y, :yerr],
        :unit_time => u"d",
        :timezero => 58541.5,
        :irregular => true
    )

    # test the assertion
    @test_throws "AssertionError: The metadata and the time stamps disagree on the time sampling. \n        metadata:true != regular_sampling false" TimeSeries(t, [x σ], metadata)

    # initialising the RegularTimeSeries with just data and times
    metadata[:irregular] = false
    ts = TimeSeries(t, [x σ], metadata)
    # check the types of the data
    @test ts.time isa RegularTimeStamps
    @test ts.data isa TimeSeriesData
    @test ts.metadata isa Dict

    # check the values in the RegularTimeSeries
    @test ts.time.Δt == 0.52
    @test ts.time.duration == t[end]
    @test ts.time.unit == u"d"
    @test ts.time.timezero == 58541.5

    # check the values in the TimeSeriesData
    @test ts.data.data == [x σ]
    @test ts.data.columns == [:y, :yerr]
    @test ts.data.units == u"erg/cm^2/s"

    # # check the values in the metadata
    @test ts.metadata[:irregular] == false
    @test ts.metadata[:unit_time] == u"d"
    @test ts.metadata[:timezero] == 58541.5
    @test ts.metadata[:columns] == [:y, :yerr]
    return @test ts.metadata[:units] == u"erg/cm^2/s"
end

function test_IrregularTimeSeries()
    t, x, σ, _, _, _, _ = setup_timeseries()

    # initialising the IrregularTimeStamps with just data and times
    @test_throws "The time stamps are not sorted in ascending order." TimeSeries(t, [x σ])
    t = sort(t)
    @test_throws "The time stamps and the data do not have the same length." TimeSeries(t, [x[2:end] σ[2:end]])

    ts = TimeSeries(t, [x σ])

    # check the types of the data
    @test ts.time isa IrregularTimeStamps
    @test ts.data isa TimeSeriesData
    @test ts.metadata isa Dict

    # check the values in the IrregularTimeStamps
    @test ts.time.Δt == minimum(diff(t))
    @test ts.time.duration == t[end] - t[1]
    @test ts.time.unit == u"s"
    @test ts.time.timezero == t[1]

    # check the values in the TimeSeriesData
    @test ts.data.data == [x σ]
    @test ts.data.columns == [:y1, :y2]
    @test ts.data.units == [NoUnits, NoUnits]

    # check the values in the metadata
    @test ts.metadata[:irregular] == true
    @test ts.metadata[:unit_time] == u"s"
    @test ts.metadata[:timezero] == t[1]
    @test ts.metadata[:columns] == [:y1, :y2]
    return @test ts.metadata[:units] == [NoUnits, NoUnits]
end

function test_IrregularTimeSeries_metadata()
    t, x, σ, _, _, _, _ = setup_timeseries()
    t = sort(t)

    metadata = Dict(
        :name => "Time series",
        :description => "A time series for testing",
        :units => u"erg/cm^2/s",
        :columns => [:y, :yerr],
        :unit_time => u"K",
        :timezero => 58541.5,
        :irregular => false
    )

    # test the assertions
    @test_throws "AssertionError: The metadata and the time stamps disagree on the time sampling. \n        metadata:false != regular_sampling true" TimeSeries(t, [x σ], metadata)
    metadata[:irregular] = true
    @test_throws "The unit of the time stamps does not have a Time (Unitful.𝐓) dimension. https://painterqubits.github.io/Unitful.jl/stable/defaultunits/#Base-dimensions" TimeSeries(t, [x σ], metadata)
    metadata[:unit_time] = u"d"

    # initialising the IrregularTimeSeries with just data and times
    metadata[:irregular] = true
    ts = TimeSeries(t, [x σ], metadata)
    # check the types of the data
    @test ts.time isa IrregularTimeStamps
    @test ts.data isa TimeSeriesData
    @test ts.metadata isa Dict

    # check the values in the IrregularTimeStamps
    @test ts.time.Δt == minimum(diff(t))
    @test ts.time.duration == t[end] - t[1]
    @test ts.time.unit == u"d"
    @test ts.time.timezero == 58541.5

    # check the values in the TimeSeriesData
    @test ts.data.data == [x σ]
    @test ts.data.columns == [:y, :yerr]
    @test ts.data.units == u"erg/cm^2/s"

    # # check the values in the metadata
    @test ts.metadata[:irregular] == true
    @test ts.metadata[:unit_time] == u"d"
    @test ts.metadata[:timezero] == 58541.5
    @test ts.metadata[:columns] == [:y, :yerr]
    return @test ts.metadata[:units] == u"erg/cm^2/s"
end

function test_fill_gaps()
    rng = MersenneTwister(123)

    t, x, σ = generate_regulartimeseries()

    x = abs.(x)
    # if nothing to interpolate or fill
    @test_throws "No missing times found check the time array" fill_gaps(rng, t, x, σ, Δt = 0.52)
    # create gaps in the time series
    index_gaps = [4, 5, 20, 29, 50, 81, 82, 82, 105, 134, 150]
    # remove the gaps
    t = t[setdiff(1:end, index_gaps)]
    x = x[setdiff(1:end, index_gaps)]
    σ = σ[setdiff(1:end, index_gaps)]

    # fill the gaps
    t_filled, x_filled, σ_filled = fill_gaps(rng, t, x, σ, Δt = 0.52)
    @test Tonari.regular_time_sampling(collect(t_filled)) == true

    return t_filled, x_filled = fill_gaps(rng, t, x, Δt = 0.52)

    # fill the gaps without randomising
    #t_filled, x_filled = fill_gaps(rng, t, x, σ, Δt = 0.52, randomise_values = false)

    # fill the gaps with a different time step
    #t_filled, x_filled, σ_filled = fill_gaps(rng, t, x, σ, randomise_values = false)
end

@testset "Time series tests" begin
    test_sanity_check()
    @testset "RegularTimeSeries" begin
        test_RegularTimeSeries()
        test_RegularTimeSeries_metadata()
    end
    @testset "IrregularTimeSeries" begin
        test_IrregularTimeSeries()
        test_IrregularTimeSeries_metadata()
    end
    @testset "Fill gaps" begin
        test_fill_gaps()
    end
end
