module IsoMix
import StaticArrays # for fast matrix inversion with the \ operator in `fractions` function
import Polyester
import Random

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

