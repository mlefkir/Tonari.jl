using Tonari
using Random
using StatsBase
using Test

@testset "Tonari.jl" begin
    include("Aqua.jl")
    include("test_simulations.jl")
    include("test_psd.jl")
    include("test_periodogram.jl")
    include("test_crossperiodogram.jl")
    include("test_timeseries.jl")
    include("test_iccf.jl")
end
