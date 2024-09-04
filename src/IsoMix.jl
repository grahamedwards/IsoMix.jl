module IsoMix
import StaticArrays, Polyester, Random

export IsoMixType, Data, DistributionDraws, Calculations
include("containers.jl")

export mix, mix!, fractions
include("mix.jl")

export loglikelihood
include("statistics.jl")

export update
include("metropolis-helpers.jl")

export mixtropolis
include("metropolis.jl")


module Examples include("examples.jl") end
end

