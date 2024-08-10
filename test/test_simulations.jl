using Test, Tonari, Random

function setup_model()
	SingleBendingPowerLaw(7e2, 0.13, 9e-2, 3.42)
end


## Test the initialisation of the Simulation object
function init_simu_regular()
	model = setup_model()
	T = 200.0
	Δt = 0.1
	sim = Simulation(model, T, Δt)
	@test sim.model == model
	@test sim.T == T - Δt
	@test sim.Δt == Δt
	@test sim.t == 0:Δt:T-Δt
end

function init_simu_regular_bis()
	model = setup_model()
	T, Δt = 302.13, .0321
	sim = Simulation(model, T, Δt)
	@test sim.model == model
	@test sim.T == T - Δt
	@test sim.Δt == Δt
	@test sim.t == 0:Δt:T-Δt
end

function init_simu_irregular()
    model = setup_model()
    rng = MersenneTwister(1234)
    t = rand(rng, 0:0.0213:236.53, 105)
    t = sort(t)
    t = unique(t)
    t .-= t[1]

    sim = Simulation(model,t)
    @test sim.model == model
    @test sim.T == t[end]
    @test sim.Δt == minimum(diff(t))
    @test sim.t == t
end

function init_simu_regular_error()
    model = setup_model()
    T, Δt = 2.13, 334.0321
    @test_throws "The sampling period Δt must be less than the duration T" Simulation(model, T, Δt)
end

function init_simu_regular_extend()
    model = setup_model()
	T, Δt = 302.13, .0321
	sim = Simulation(model, T, Δt,15.5,14.4)
	@test sim.model == model
	@test sim.T == T - Δt
	@test sim.Δt == Δt
	@test sim.t == 0:Δt:T-Δt
    @test sim.S_high == 15.5
    @test sim.S_low == 14.4
end

function init_simu_irregular_unsorted()
    model = setup_model()
    rng = MersenneTwister(1234)
    t = rand(rng, 0:0.0213:236.53, 105)
    t = unique(t)
    t .-= t[1]

    @test_throws "The time vector must be sorted" Simulation(model,t)
end


function init_simu_irregular_extend()
    model = setup_model()
    rng = MersenneTwister(1234)
    t = rand(rng, 0:0.0213:236.53, 105)
    t = sort(t)
    t = unique(t)
    t .-= t[1]

    sim = Simulation(model,t,15.5,14.4)
    @test sim.model == model
    @test sim.T == t[end]
    @test sim.Δt == minimum(diff(t))
    @test sim.t == t
    @test sim.S_high == 15.5
    @test sim.S_low == 14.4
end


@testset "Simulations" begin
	@testset "Initialisation" begin
        @testset "Regular" begin
            init_simu_regular()
            init_simu_regular_bis()
            init_simu_regular_error()
            init_simu_regular_extend()
        end
        @testset "Irregular" begin
            init_simu_irregular()
            init_simu_irregular_extend()
            init_simu_irregular_unsorted()
        end
    end
end
