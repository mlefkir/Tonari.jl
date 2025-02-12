using Test, Tonari, Random

function test_init_cross_spectrum_timelag()
    p1 = SingleBendingPowerLaw(1.0, 0.63, 10.0e-2, 3.22)
    p2 = SingleBendingPowerLaw(1.0, 0.78, 13.0e-2, 3.81)
    Δϕ = ConstantTimeLag(100.35)
    cp = CrossSpectralDensity(p1, p2, Δϕ)
    @test cp.𝓟₁ == p1
    @test cp.𝓟₂ == p2
    return @test cp.Δφ == Δϕ
end

function test_init_cross_spectrum_phaselag()
    p1 = SingleBendingPowerLaw(1.0, 0.63, 10.0e-2, 3.22)
    p2 = SingleBendingPowerLaw(1.0, 0.78, 13.0e-2, 3.81)
    Δϕ = ConstantPhaseLag(35.35, 2.3)
    cp = CrossSpectralDensity(p1, p2, Δϕ)
    @test cp.𝓟₁ == p1
    @test cp.𝓟₂ == p2
    return @test cp.Δφ == Δϕ
end

function test_cross_periodogram()
    p1 = SingleBendingPowerLaw(1.0, 0.63, 10.0e-2, 3.22)
    p2 = SingleBendingPowerLaw(1.0, 0.78, 13.0e-2, 3.81)
    Δϕ = ConstantTimeLag(100.35)
    model = CrossSpectralDensity(p1, p2, Δϕ)
    T = 200.0
    Δt = 0.1
    sim = Simulation(model, T, Δt)
    n = 20
    t, x, σ = sample(MersenneTwister(1234), sim, n)
    f, γ², γ²_corrected, Δφ, γ²_err, γ²_corrected_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, N₁, N₂, n = cross_periodogram(t, x[1], x[2], σ[1], σ[2])
    @test length(f) == length(γ²) == length(γ²_corrected) == length(Δφ) == length(γ²_err) == length(γ²_corrected_err) == length(Δφ_err) == length(Δτ) == length(Δτ_err) == length(P̄₁) == length(P̄₂) == length(n)
    f, C̄, P̄₁, P̄₂, C = cross_periodogram(t, x[1], x[2], σ[1], σ[2], compute_coherence = false)
    @test length(f) == length(C̄) == length(P̄₁) == length(P̄₂)
    f, γ², γ²_corrected, Δφ, γ²_err, γ²_corrected_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, N₁, N₂, n = cross_periodogram(t, x[1], x[2], σ[1], σ[2], apply_end_matching = true, subtract_mean = false)
    return f, γ², γ²_corrected, Δφ, γ²_err, γ²_corrected_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, N₁, N₂, n = cross_periodogram(t, x[1], x[2], σ[1], σ[2], subtract_mean = false)


end

function test_cross_periodogram_noerr()
    p1 = SingleBendingPowerLaw(1.0, 0.63, 10.0e-2, 3.22)
    p2 = SingleBendingPowerLaw(1.0, 0.78, 13.0e-2, 3.81)
    Δϕ = ConstantTimeLag(100.35)
    model = CrossSpectralDensity(p1, p2, Δϕ)
    T = 200.0
    Δt = 0.1
    sim = Simulation(model, T, Δt)
    n = 20
    t, x, _ = sample(MersenneTwister(1234), sim, n)
    f, γ², Δφ, γ²_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, C̄ = cross_periodogram(t, x[1], x[2])
    return @test length(f) == length(γ²) == length(Δφ) == length(γ²_err) == length(Δφ_err) == length(Δτ) == length(Δτ_err) == length(P̄₁) == length(P̄₂) == length(C̄)

end

function test_cross_periodogram_can_fail()
    p1 = SingleBendingPowerLaw(1.0, 0.63, 10.0e-2, 3.22)
    p2 = SingleBendingPowerLaw(1.0, 0.78, 13.0e-2, 3.81)
    Δϕ = ConstantTimeLag(100.35)
    model = CrossSpectralDensity(p1, p2, Δϕ)
    T = 200.0
    Δt = 0.1
    sim = Simulation(model, T, Δt)
    n = 1
    t, x, _ = sample(MersenneTwister(1234), sim, n)
    return @test_throws "Cannot compute coherence or lags without several segments of the full time series." f, γ², Δφ, γ²_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, C̄ = cross_periodogram(t, x[1], x[2])

end


@testset "Modelling" begin
    test_init_cross_spectrum_timelag()
    test_init_cross_spectrum_phaselag()
end

@testset "cross_periodogram" begin
    test_cross_periodogram()
    test_cross_periodogram_noerr()
end
