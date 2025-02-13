using Test, Tonari, Random, Distributions

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

@testset "Power spectrum model" begin
    test_init_model_sum()
end
