using Test, Tonari, Random

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

    @test_logs (:info, "Sorting the time series")  time_series_sanity_checks(t₂, x₂)
    @test_throws "The time series contains infinities or NaNs" time_series_sanity_checks(t₁_sort, x₁)
    return @test_throws "Uncertainty on the time series contains infinities or NaNs" time_series_sanity_checks(sort(t₂), x₂, σ₂)
end


@testset "Time series tests" begin
    test_sanity_check()
end
