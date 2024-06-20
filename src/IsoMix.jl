module IsoMix
import StaticArrays, Polyester, Random

export IsoMixType
export Component, Component1, Component2, Component3
export Datum, Norm, logNorm, Unf, Constant, Unconstrained
export Data, Data1, Data2, Data3
export Prior, Prior2, Prior3
export System, System2, System3 
export Fraction, Fraction2, Fraction3
export Model, Model1, Model2, Model3
export Measurements, Measurements1, Measurements2, Measurements3
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

