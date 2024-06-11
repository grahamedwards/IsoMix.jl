module IsoMix

export IsoMixType
export Component, countcomponents
export Datum, Norm, logNorm, Unf, Constant, Data, Data1, Data2, Data3, Prior
export System, System2, System3 
export Fraction, Fraction2, Fraction3
export Model, Model1, Model2, Model3
include("containers.jl")


export mix, mix!
include("mix.jl")

end
