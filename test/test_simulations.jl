using Test, Tonari, Random, Distributions, DelimitedFiles

function setup_model()
	SingleBendingPowerLaw(7e2, 0.13, 9e-2, 3.42)
end

function setup_model_multivariate_time()
	p1 = SingleBendingPowerLaw(1.0, 0.63, 10e-2, 3.22)
	p2 = SingleBendingPowerLaw(1.0, 0.78, 13e-2, 3.81)
	Δϕ = ConstantTimeLag(100.35)
	cp = CrossSpectralDensity(p1, p2, Δϕ)
end

function setup_model_multivariate_phase()
	p1 = SingleBendingPowerLaw(1.0, 0.63, 10e-2, 3.22)
	p2 = SingleBendingPowerLaw(1.0, 0.78, 13e-2, 3.81)
	Δϕ = ConstantPhaseLag(35.35, 2.3)
	cp = CrossSpectralDensity(p1, p2, Δϕ)
end


function setup_model_multivariate_same_psd()
	p1 = SingleBendingPowerLaw(1.0, 0.63, 10e-2, 3.22)
	Δϕ = ConstantPhaseLag(35.35, 2.3)
	cp = CrossSpectralDensity(p1, p1, Δϕ)
end

## Test the initialisation of the Simulation object
function init_simu_regular()
	model = setup_model()
	T = 200.0
	Δt = 0.1
	sim = Simulation(model, T, Δt)
	@test sim.model == model
	@test sim.T == T
	@test sim.Δt == Δt
	@test sim.t == 0:Δt:T
end

function init_simu_regular_bis()
	model = setup_model()
	T, Δt = 302.13, 0.0321
	sim = Simulation(model, T, Δt)
	@test sim.model == model
	@test sim.T == T
	@test sim.Δt == Δt
	@test sim.t == 0:Δt:T
end

function init_simu_irregular()
	model = setup_model()
	rng = MersenneTwister(1234)
	t = rand(rng, 0:0.0213:236.53, 105)
	t = sort(t)
	t = unique(t)
	t .-= t[1]

	sim = Simulation(model, t)
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
	T, Δt = 302.13, 0.0321
	sim = Simulation(model, T, Δt, 15.5, 14.4)
	@test sim.model == model
	@test sim.T == T
	@test sim.Δt == Δt
	@test sim.t == 0:Δt:T
	@test sim.S_high == 15.5
	@test sim.S_low == 14.4
end

function init_simu_irregular_unsorted()
	model = setup_model()
	rng = MersenneTwister(1234)
	t = rand(rng, 0:0.0213:236.53, 105)
	t = unique(t)
	t .-= t[1]

	@test_throws "The time vector must be sorted" Simulation(model, t)
end


function init_simu_irregular_extend()
	model = setup_model()
	rng = MersenneTwister(1234)
	t = rand(rng, 0:0.0213:236.53, 105)
	t = sort(t)
	t = unique(t)
	t .-= t[1]

	sim = Simulation(model, t, 15.5, 14.4)
	@test sim.model == model
	@test sim.T == t[end]
	@test sim.Δt == minimum(diff(t))
	@test sim.t == t
	@test sim.S_high == 15.5
	@test sim.S_low == 14.4
end

function init_simu_irregular_load()
	model = setup_model()
	rng = MersenneTwister(1234)
	t = vec(readdlm("data/sampling_1yr_f9.txt"))
	t = sort(t)
	t = unique(t)
	t .-= t[1]

	sim = Simulation(model, t, 15.5, 14.4)
	@test sim.model == model
	@test sim.T == t[end]
	@test sim.Δt == minimum(diff(t))
	@test sim.t == t
	@test sim.S_high == 15.5
	@test sim.S_low == 14.4
end



## Test the sampling of the Simulation object

function setUp_simu_regular()
	model = setup_model()
	T = 200.0
	Δt = 0.1
	sim = Simulation(model, T, Δt)
	return sim, T, Δt
end
function sample_simu_regular()

	sim, T, Δt = setUp_simu_regular()
	rng = MersenneTwister(306)
	tx, x, σ = sample(rng, sim)
	N = convert(Int, T / Δt)
	@test tx ≈ sim.t[1:end-1]
	@test size(x) == size(σ)
	@test size(x) == (N, 1)
end

function sample_simu_regular_n()
	sim, T, Δt = setUp_simu_regular()
	n = 3
	rng = MersenneTwister(124)
	tx, x, σ = sample(rng, sim, n)
	N = convert(Int, T / Δt)
	@test tx ≈ sim.t[1:end-1]
	@test size(x) == size(σ)
	@test size(x) == (N, n)
end

function sample_simu_regular_n_nosplit()
	sim, T, Δt = setUp_simu_regular()
	n = 3
	rng = MersenneTwister(216)
	tx, x, σ = sample(rng, sim, n, split_long = false)
	N = convert(Int, T / Δt)
	@test tx ≈ sim.t[1:end-1]
	@test size(x) == size(σ)
	@test size(x) == (N, n)
end

function sample_simu_regular_poisson()
	sim, T, Δt = setUp_simu_regular()
	n = 3
	rng = MersenneTwister(816)
	tx, x, σ = sample(rng, sim, n, poisson = true)
	N = convert(Int, T / Δt)
	@test tx ≈ sim.t[1:end-1]
	@test size(x) == size(σ)
	@test size(x) == (N, n)
	@test all(x .>= 0)
end

function sample_simu_regular_expo()
	sim, T, Δt = setUp_simu_regular()
	n = 3
	rng = MersenneTwister(936)
	tx, x, σ = sample(rng, sim, n, exponentiate = true)
	N = convert(Int, T / Δt)
	@test tx ≈ sim.t[1:end-1]
	@test size(x) == size(σ)
	@test size(x) == (N, n)
	@test all(x .>= 0)
end

function sample_simu_regular_error()
	sim, T, Δt = setUp_simu_regular()
	n = 3
	rng = MersenneTwister(2346)
	N = convert(Int, T / Δt)
	σ = rand(rng, N)
	tx, x, σs = sample(rng, sim, n, σₓ = σ)
	@test size(x, 1) == size(σ, 1)
	@test σ == σs
end

function sample_simu_irregular_load()
	model = setup_model()
	rng = MersenneTwister(1234)
	t = vec(readdlm("data/sampling_1yr_f9.txt"))
	t = sort(t)
	t = unique(t)
	t .-= t[1]

	sim = Simulation(model, t, 15.5, 14.4)
	n = 3
	tx, x, σs = sample(rng, sim, n)
	@test size(x, 1) == size(σs, 1)
	# @test size(x, 2) == (length(t)-1, n)
end

function sample_simu_irregular_extend()
	model = setup_model()
	rng = MersenneTwister(1234)
	t = rand(rng, 0:0.0213:236.53, 105)
	t = sort(t)
	t = unique(t)
	t .-= t[1]
	n = 3
	sim = Simulation(model, t, 15.5, 14.4)
	tx, x, σs = sample(rng, sim, n)
end


# Test the sampling of the Simulation object with a cross-spectral density
function sample_simu_cross_spectral_time()
	model = setup_model_multivariate_time()
	T = 200.0
	Δt = 0.1
	sim = Simulation(model, T, Δt)
	n = 4
	tx, x, σ = sample(MersenneTwister(1234), sim, n)

	@test size(x,1) == 2
	@test size(σ,1) == 2
	@test size(x[1]) == size(x[2])
	@test size(x[1]) == size(σ[1])
	@test size(x[1]) == size(σ[2])
	@test size(x[1],1) == size(tx,1)
end

function sample_simu_cross_spectral_phase()
	model = setup_model_multivariate_phase()
	T = 200.0
	Δt = 0.1
	sim = Simulation(model, T, Δt)
	n = 4
	tx, x, σ = sample(MersenneTwister(1234), sim, n)

	@test size(x,1) == 2
	@test size(σ,1) == 2
	@test size(x[1]) == size(x[2])
	@test size(x[1]) == size(σ[1])
	@test size(x[1]) == size(σ[2])
	@test size(x[1],1) == size(tx,1)
end
function sample_simu_cross_spectral_same_psd()
	model = setup_model_multivariate_same_psd()
	T = 200.0
	Δt = 0.1
	sim = Simulation(model, T, Δt)
	n = 4
	t, x, σ = sample(MersenneTwister(1234), sim, n)

	@test size(x,1) == 2
	@test size(σ,1) == 2

	@test size(x[1]) == size(x[2])
	@test size(x[1]) == size(σ[1])
	@test size(x[1]) == size(σ[2])
	@test size(x[1],1) == size(t,1)
end

function sample_simu_cross_spectral_no_rand()
	model = setup_model_multivariate_same_psd()
	T = 200.0
	Δt = 0.1
	sim = Simulation(model, T, Δt)
	n = 4
	t, x = sample(MersenneTwister(1234), sim, n,randomise_values = false)

	@test size(x[1],1) == size(t,1)
	@test size(x,1) == 2
	@test size(x[1]) == size(x[2])
end



@testset "Simulations" begin
	@testset "Univariate" begin
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
				init_simu_irregular_load()
			end
		end
		@testset "Sampling" begin
			@testset "Regular" begin
				sample_simu_regular()
				sample_simu_regular_n()
				sample_simu_regular_n_nosplit()
				sample_simu_regular_poisson()
				sample_simu_regular_expo()
				sample_simu_regular_error()
			end
			@testset "Irregular" begin
				sample_simu_irregular_load()
				sample_simu_irregular_extend()
			end
		end
	end
	@testset "Multivariate" begin
		@testset "Initialisation" begin
			@testset "Regular" begin
				sample_simu_cross_spectral_time()
				sample_simu_cross_spectral_phase()
				sample_simu_cross_spectral_same_psd()
				sample_simu_cross_spectral_no_rand()
			end
		end
	end
end
