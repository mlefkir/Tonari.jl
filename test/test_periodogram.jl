using Test, Tonari, Random

function setUp(n = 1, N = 2000)
    model = SingleBendingPowerLaw(7.0e2, 0.13, 9.0e-2, 3.42)
    Δt = 0.1
    T = N * Δt
    rng = MersenneTwister(12)
    sim = Simulation(model, T, Δt)
    t, x, σ = sample(rng, sim, n)
    return t, x, σ
end

function test_periodogram_single()
    t, x, σ = setUp()
    f, I = periodogram(t, x)
    return @test length(f) == length(I)
end

function test_periodogram_average()
    t, x, σ = setUp(100)
    f, I = periodogram(t, x)
    return @test length(f) == length(I)
end

function test_periodogram_average_endmatch()
    t, x, σ = setUp(100)
    f, I = periodogram(t, x, apply_end_matching = true)
    return @test length(f) == length(I)
end

function test_periodogram_single_endmatch()
    t, x, σ = setUp()
    f, I = periodogram(t, x, apply_end_matching = true)
    return @test length(f) == length(I)
end


function test_periodogram_nonequi()
    model = SingleBendingPowerLaw(7.0e2, 0.13, 9.0e-2, 3.42)
    rng = MersenneTwister(12)
    t = rand(rng, 0:0.0213:236.53, 1205)
    t = sort(t)
    t = unique(t)
    sim = Simulation(model, t)
    t, x, σ = sample(rng, sim)
    return @test_throws "The time stamps are not equally spaced! Check your time array!"  periodogram(t, x)
end

@testset "Periodogram estimation" begin
    test_periodogram_single()
    test_periodogram_single_endmatch()
    test_periodogram_average()
    test_periodogram_average_endmatch()
    test_periodogram_nonequi()
end
