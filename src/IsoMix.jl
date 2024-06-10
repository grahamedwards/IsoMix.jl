module IsoMix

export Component, countcomponents
export System, System2, System3 
export Prior, Norm, logNorm, Unf, Constant, Data, Data2, Data3
export Fraction, Fraction2, Fraction3
export Model, Model2, Model3
include("containers.jl")


export mix, mix!
include("mix.jl")

end
