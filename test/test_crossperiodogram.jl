using Test, Tonari, Random

function test_init_cross_spectrum_timelag()
	p1 = SingleBendingPowerLaw(1.0, 0.63, 10e-2, 3.22)
	p2 = SingleBendingPowerLaw(1.0, 0.78, 13e-2, 3.81)
	Δϕ = ConstantTimeLag(100.35)
	cp = CrossPeriodogram(p1, p2, Δϕ)
	@test cp.𝓟₁ == p1
	@test cp.𝓟₂ == p2
	@test cp.Δφ == Δϕ
end

function test_init_cross_spectrum_phaselag()
	p1 = SingleBendingPowerLaw(1.0, 0.63, 10e-2, 3.22)
	p2 = SingleBendingPowerLaw(1.0, 0.78, 13e-2, 3.81)
	Δϕ = ConstantPhaseLag(.35)
	cp = CrossPeriodogram(p1, p2, Δϕ)
	@test cp.𝓟₁ == p1
	@test cp.𝓟₂ == p2
	@test cp.Δφ == Δϕ
end


@testset "CrossPeriodogram" begin
    test_init_cross_spectrum_timelag()
    test_init_cross_spectrum_phaselag()
end