using Test, Tonari, Random

function generate_lagged_ts()
    p1 = SingleBendingPowerLaw(1.0, 0.63, 10.0e-2, 3.22)
    p2 = SingleBendingPowerLaw(1.0, 0.78, 13.0e-2, 4.81)
    Δϕ = ConstantTimeLag(18.57)
    cs = CrossSpectralDensity(p1, p2, Δϕ)

    rng = MersenneTwister(0)
    T, Δt = 504.2, 0.5
    simu = Simulation(cs, T, Δt)
    t, y, yerr = sample(rng, simu, 1, error_size = 0.15)
    N1, N2 = 511, 680
    # sample indices
    p1 = sort(sample(rng, 1:length(t), N1, replace = false))
    p2 = sort(sample(rng, 1:length(t), N2, replace = false))

    t₁, y₁, σ₁ = t[p1], y[1][p1, 1], yerr[1][p1, 1]
    t₂, y₂, σ₂ = t[p2], y[2][p2, 1], yerr[2][p2, 1]
    return t₁, y₁, t₂, y₂, σ₁, σ₂
end

""" test the iccf function """
function test_iccf()
    t₁, y₁, t₂, y₂, σ₁, σ₂ = generate_lagged_ts()


    τ_list, r = cross_correlate(t₁, y₁, t₂, y₂, σ₁ = σ₁, σ₂ = σ₂)
    @test length(r) == length(τ_list)
    τ_in = -40:0.5:40
    τ_list2, r = cross_correlate(t₁, y₁, t₂, y₂, τ_list = τ_in, σ₁ = σ₁, σ₂ = σ₂)
    @test length(r) == length(τ_in)

    rng = MersenneTwister(0)
    #τ_in = -40:0.5:40
    τ_list, r, q = cross_correlate(rng, t₁, y₁, t₂, y₂, σ₁ = σ₁, σ₂ = σ₂, compute_errors = true, n_simulations = 450)
    @test length(q) == 450
    return @test length(r) == length(τ_list)
end

@testset "ICCF" begin
    test_iccf()
end
