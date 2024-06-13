using IsoMix
using Test
using StableRNGs

include("silence.jl")


## Some helpful functions
μ(x) = sum(x)/length(x)
Base.isapprox(x::T,y::T) where T <: AbstractArray = prod( x .≈ y) # Returns true if x[i] ≈ y[i] for all i.

@testset "Containers and constructors" begin include("containers.jl") end
@testset "Mixing math" begin include("mix.jl") end
