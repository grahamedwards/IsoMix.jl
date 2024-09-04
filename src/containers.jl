"""
    IsoMixType

Supertype containing all custom types and structs in the IsoMix package:

Direct subtypes: [`EmDraw`](@ref), [`Endmember`](@ref), [`Measurements`](@ref), [`Datum`](@ref), [`Fraction`](@ref), [`Measurements`](@ref), [`Model`](@ref),[`Prior`](@ref), [`SysDraw`](@ref) 

---
```
IsoMixType
├─ Calculations
│  ├─ Fraction
│  │  ├─ Fraction2
│  │  └─ Fraction3
│  └─ Model
│     ├─ Model1
│     ├─ Model2
│     └─ Model3
├─ Data
│  ├─ Datum
│  │  ├─ Constant
│  │  ├─ Norm
│  │  ├─ Unconstrained
│  │  ├─ Unf
│  │  └─ logNorm
│  ├─ Endmember
│  │  ├─ Endmember1
│  │  ├─ Endmember2
│  │  └─ Endmember3
│  ├─ Measurements
│  │  ├─ Measurements1
│  │  ├─ Measurements2
│  │  └─ Measurements3
│  └─ Prior
│     ├─ Prior2
│     └─ Prior3
└─ DistributionDraws
   ├─ EmDraw
   │  ├─ EmDraw1
   │  ├─ EmDraw2
   │  └─ EmDraw3
   └─ SysDraw
      ├─ SysDraw2
      └─ SysDraw3
```

---
Generated with 

    using AbstractTrees,IsoMix; AbstractTrees.children(d) = subtypes(d); print_tree(IsoMixType)

"""
abstract type IsoMixType end

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""

    Data <: IsoMixType

Abstract supertype for types that contain measured data, represented by `Datum` instances. `Prior` instances contain 2-3 `Endmember` instances, and both `Endmember or `Measurements` instances contain `Datum` instances.

see also: [`Datum`](@ref), [`Endmember`](@ref), [`Prior`](@ref), [`Measurements`](@ref)

"""
abstract type Data <: IsoMixType end

include("containers/data.jl")
export Datum, Norm, logNorm, Unf, Constant, Unconstrained
export Endmember, Endmember1, Endmember2, Endmember3
export Measurements, Measurements1, Measurements2, Measurements3
export Prior, Prior2, Prior3

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



"""

    DistributionDraws <: IsoMixType

Abstract supertype for types that contain values drawn from [`Data`](@ref) distributions of `Endmember` and `Prior` instances. `SysDraw` instances contain 2-3 `EmDraw` instances, which contain discrete draws from a corresponding `Endmember` instance.

see also: [`EmDraw`](@ref), [`SysDraw`](@ref)

"""
abstract type DistributionDraws <: IsoMixType end 

include("containers/distribution-draws.jl")
export EmDraw, EmDraw1, EmDraw2, EmDraw3
export SysDraw, SysDraw2, SysDraw3 



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



"""

    Calculations <: IsoMixType

Abstract supertype for types used specifically in chemical and isotopic model simulations. 

see also: [`Fraction`](@ref), [`Model`](@ref)

"""
abstract type Calculations<: IsoMixType end

include("containers/calculations.jl")
export Fraction, Fraction2, Fraction3
export Model, Model1, Model2, Model3



"""

    IsoMix.countcomponents(c::EmDraw)

Count the number of elemental dimensions (1-3) of `c`.

# Example

    julia> IsoMix.countcomponents(EmDraw(1,2,3,4,5,6))
    3

"""
countcomponents(::C) where C <: EmDraw = div(fieldcount(C),2)