using Tonari
using Random
using StatsBase
using Test

@testset "Tonari.jl" begin
    include("test_simulations.jl")
    include("test_psd.jl")
    include("test_periodogram.jl")
    include("test_crossperiodogram.jl")
end
