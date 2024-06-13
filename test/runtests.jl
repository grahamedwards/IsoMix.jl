using IsoMix
using Test
using StableRNGs

include("silence.jl")

μ(x) = sum(x)/length(x)

@testset "Containers and constructors" begin include("containers.jl") end
@testset "Mixing math" begin include("mix.jl") end
