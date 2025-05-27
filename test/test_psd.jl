using Test, Tonari, Random, Distributions

function test_SingleBendingPowerLaw()
    PS = SingleBendingPowerLaw(0.3, 0.02, 2.93)
    f = 10 .^ range(-3, stop = 2, length = 1000)
    return @test PS(f) == (f / 0.02) .^ (-0.3) ./ (1 .+ (f / 0.02) .^ (2.93 - 0.3))
end

function test_DoubleBendingPowerLaw()
    DS = DoubleBendingPowerLaw(0.3, 0.02, 1.4, 10.2, 2.93)
    f = 10 .^ range(-3, stop = 3, length = 1000)
    return @test DS(f) == (f / 0.02) .^ (-0.3) ./ (1 .+ (f / 0.02) .^ (1.4 - 0.3)) ./ (1 .+ (f / 10.2) .^ (2.93 - 1.4))
end


function test_QPO()
    Q = QPO(1.4, 0.2, 6.3)
    f = 10 .^ range(-3, stop = 3, length = 1000)
    return @test Q(f) == 1.4 * 0.2^4 ./ ((f .^ 2 .- 0.2^2) .^ 2 .+ f .^ 2 .* 0.2^2 / 6.3^2)
end

function test_PowerLaw()
    PL = PowerLaw(1.4)
    f = 10 .^ range(-3, stop = 3, length = 1000)
    return @test PL(f) == f .^ (-1.4)
end

function setUp_models()
    m1 = SingleBendingPowerLaw(7.0e2, 0.13, 9.0e-2, 3.42)
    m2 = DoubleBendingPowerLaw(0.13, 9.0e-2, 3.42, 0.5, 0.1)
    m3 = Lorentzian(7.0e2, 0.13, 9.0e-2)
    return m1, m2, m3
end

function test_init_model_sum()
    m1, m2, m3 = setUp_models()
    m = m1 + m2 + m3

    # test the values of the model
    @test m.psd == [m1, m2, m3]
    @test ((m1 + m2) + m3).psd == [m1, m2, m3]
    @test (m1 + (m2 + m3)).psd == [m1, m2, m3]
    @test (m3 + (m1 + m2)).psd == [m3, m1, m2]
    f = 10 .^ range(-3, 2, length = 1000)
    # test the evaluation of the model
    return @test m(f) == m1(f) + m2(f) + m3(f)
end

function test_separate_psd_models()
    m1, m2 = SingleBendingPowerLaw(7.0e2, 0.13, 9.0e-2, 3.42), QPO(1.03, 0.2, 9.4)

    # check that if the model is a continuum model it returns the model
    u1, u2 = separate_psd(m1)
    @test u1 == m1
    @test isnothing(u2)

    # if we have a QPO only returns QPO
    u1, u2 = separate_psd(m2)
    @test u2 == m2
    @test isnothing(u1)

    # if we have a sum of continuum and QPO
    u1, u2 = separate_psd(m1 + m2)
    @test u1 == m1
    @test u2[1] == m2

    # if we have two Lorentzian
    u1, u2 = separate_psd(m2 + m2)
    @test isnothing(u1)
    @test u2 isa Vector
    return @test u2[1] == m2 && u2[2] == m2

    # if we have BendingPowerLaws
    u1, u2 = separate_psd(m1 + m1)
    @test isnothing(u2)
    @test u1 isa SumOfPowerSpectralDensity
    @test u1.psd[1] == u1.psd[2] == m1
end

@testset "Power spectrum model" begin
    test_DoubleBendingPowerLaw()
    test_SingleBendingPowerLaw()
    test_QPO()
    test_PowerLaw()
    test_init_model_sum()
    test_separate_psd_models()
end
