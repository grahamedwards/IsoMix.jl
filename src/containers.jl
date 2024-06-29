"""
    IsoMixType

Supertype containing all custom types and structs in the IsoMix package:

Direct subtypes: [`EmDraw`](@ref), [`Endmember`](@ref), [`Measurements`](@ref), [`Datum`](@ref), [`Fraction`](@ref), [`Measurements`](@ref), [`Model`](@ref),[`Prior`](@ref), [`SysDraw`](@ref) 

---
```
IsoMixType
├─ EmDraw
│  ├─ EmDraw1
│  ├─ EmDraw2
│  └─ EmDraw3
├─ Endmember
│  ├─ Endmember1
│  ├─ Endmember2
│  └─ Endmember3
├─ Datum
│  ├─ Constant
│  ├─ Norm
│  ├─ Unf
│  └─ logNorm
├─ Fraction
│  ├─ Fraction2
│  └─ Fraction3
├─ Measurements
├─ Model
│  ├─ Model1
│  ├─ Model2
│  └─ Model3
├─ Prior
│  ├─ Prior2
│  └─ Prior3
└─ SysDraw
   ├─ SysDraw2
   └─ SysDraw3
```

---
Generated with 

    using AbstractTrees,IsoMix; AbstractTrees.children(d) = subtypes(d); print_tree(IsoMixType)

"""
abstract type IsoMixType end


"""

    DistributionDraws <: IsoMixType

Abstract supertype for types that contain values drawn from [`Data`](@ref) distributions of `Endmember` and `Prior` instances. `SysDraw` instances contain 2-3 `EmDraw` instances, which contain discrete draws from a corresponding `Endmember` instance.

see also: [`EmDraw`](@ref), [`SysDraw`](@ref)

"""
abstract type DistributionDraws <: IsoMixType end 

"""

    Data <: IsoMixType

Abstract supertype for types that contain measured data, represented by `Datum` instances. `Prior` instances contain 2-3 `Endmember` instances, and both `Endmember or `Measurements` instances contain `Datum` instances.

see also: [`Datum`](@ref), [`Endmember`](@ref), [`Prior`](@ref), [`Measurements`](@ref)

"""
abstract type Data <: IsoMixType end

## EmDraw

"""
    EmDraw <: DistributionDraws

Abstract supertype for `EmDraw_` instances, where the `_` indicates dimensionality. Represents the composition of an endmember within a natural system. 

see also: [`SysDraw`](@ref), [`EmDraw1`](@ref), [`EmDraw2`](@ref), [`EmDraw3`](@ref)

---

`EmDraw(x, cx)` -> `EmDraw1`

`EmDraw(x, cx, y, cy)` -> `EmDraw2`

`EmDraw(x, cx, y, cy, z, cz)` -> `EmDraw3`

The constructor function accepts values as a list of arguments or as keyword assignments, e.g.: `EmDraw(x=1, cx=2)`

"""
abstract type EmDraw <: DistributionDraws end
EmDraw(x::Number, cx::Number) = EmDraw1(x, cx)

EmDraw(x::Number, cx::Number, y::Number, cy::Number) = EmDraw2(x, cx, y, cy)

EmDraw(x::Number, cx::Number, y::Number, cy::Number, z::Number, cz::Number) = EmDraw3(x, cx, y, cy, z, cz)

function EmDraw(;x::Number=NaN, cx::Number=NaN, y::Number=NaN, cy::Number=NaN, z::Number=NaN, cz::Number=NaN) 
    out = if (z == z) & (cz == cz)
        EmDraw3(x, cx, y, cy, z, cz)
    elseif (y==y) & (cy==cy)
        EmDraw2(x, cx, y, cy)
    elseif (x==x) & (cx==cx)
        EmDraw1(x,cx)
    else
        error("Improper EmDraw definition")
    end
    out
end

"""
    EmDraw1 <: EmDraw

1-dimensional `EmDraw` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x

see also: [`EmDraw`](@ref)

"""
struct EmDraw1 <: EmDraw
    x::Float64
    cx::Float64
end 

"""
    EmDraw2 <: EmDraw

2-dimensional `EmDraw` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y

see also: [`EmDraw`](@ref)

"""
struct EmDraw2 <: EmDraw
    x::Float64
    cx::Float64
    y::Float64
    cy::Float64
end 

"""
    EmDraw3 <: EmDraw

3-dimensional `EmDraw` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y
`z` | Isotopic composition of z
`cz`| Concentration of z

see also: [`EmDraw`](@ref)

"""
struct EmDraw3 <: EmDraw
    x::Float64
    cx::Float64
    y::Float64
    cy::Float64
    z::Float64
    cz::Float64
end



"""

    IsoMix.countcomponents(c::EmDraw)

Count the number of elemental dimensions (1-3) of `c`.

# Example

    julia> IsoMix.countcomponents(EmDraw(1,2,3,4,5,6))
    3

"""
countcomponents(::C) where C <: EmDraw = div(fieldcount(C),2)

## SysDraw

"""

    SysDraw <: Data

Abstract supertype for `SysDraw_` instances, where the `_` indicates dimensionality. Represents a model system composition charaterized by its `EmDraw` fields.

see also: [`EmDraw`](@ref) [`SysDraw2`](@ref), [`SysDraw3`](@ref)

---
The constructor function accepts values as a list of arguments:

    SysDraw(A::EmDraw, B::EmDraw) -> SysDraw2

    SysDraw(A::EmDraw, B::EmDraw, C::EmDraw) -> SysDraw3

... or as keyword assignments, e.g.: 
    SysDraw(A=..., B=...)

"""
abstract type SysDraw{T<:EmDraw} <: DistributionDraws end
SysDraw(A::T, B::T) where T <: EmDraw = SysDraw2(A,B)
SysDraw(A::T, B::T, C::T) where T <: EmDraw = SysDraw3(A,B,C)

"""
    SysDraw2 <: SysDraw

2-dimensional `SysDraw` instance.

|Fields ||
|:--|:--|
`A` | Endmember/component A
`B` | Endmember/component B

see also: [`SysDraw`](@ref)
"""
@kwdef struct SysDraw2{T<:EmDraw} <: SysDraw{T}
    A::T
    B::T
end

"""
    SysDraw3 <: SysDraw

3-dimensional `SysDraw` instance.

|Fields ||
|:--|:--|
`A` | Endmember/component A
`B` | Endmember/component B
`C` | Endmember/component C

see also: [`SysDraw`](@ref)
"""
@kwdef struct SysDraw3{T<:EmDraw} <: SysDraw{T}
    A::T
    B::T
    C::T
end

## Datum

"""
    Datum <: Data

Abstract supertype for `Datum` instances.

see also: [`Constant`](@ref), [`Unf`](@ref), [`Norm`](@ref), [`logNorm`](@ref), [`Unconstrained`](@ref)

"""
abstract type Datum <: Data end 

"""

    Norm(m, s) <: Datum

Normally distributed [`Datum`](@ref) with mean `m` and standard deviation `s`.

"""
struct Norm <: Datum
    m::Float64
    s::Float64
end

"""

    logNorm(lm, ls) <: Datum

Log-normally distributed [`Datum`](@ref) with log-mean `lm` and log-space standard deviation `ls`.

"""
struct logNorm <: Datum
    lm::Float64
    ls::Float64
end

"""

    Unf(a, b) <: Datum

Uniformly distributed [`Datum`](@ref) with minimum `a` and maximum `b`, inclusive. I.e., for any \$x\$ in \$Unf(a, b)\$, \$x \\in [a,b]\$.

"""
struct Unf <: Datum
    a::Float64
    b::Float64
end

"""

    Constant(x) <: Datum

A [`Datum`](@ref) instance with a discrete value `x`.

"""
struct Constant <: Datum
    x::Float64
end

"""

    Unconstrained <: Datum

A [`Datum`](@ref) instance for an unconstrained variable, functionally similar to `Unf(-∞,∞)` or an unmeasured `Measurements` datum. Typically, using this value is inadvisable.

"""
struct Unconstrained <: Datum end

## Endmember

"""

    Endmember <: Data

Abstract supertype for `Endmember_` instances, where the `_` indicates dimensionality.

see also: [`Endmember1`](@ref), [`Endmember2`](@ref), [`Endmember3`](@ref)

"""
abstract type Endmember <: Data end

Endmember(x::Datum, cx::Datum) = Endmember1(x, cx)

Endmember(x::Datum, cx::Datum, y::Datum, cy::Datum) = Endmember2(x, cx, y, cy)

Endmember(x::Datum, cx::Datum, y::Datum, cy::Datum, z::Datum, cz::Datum) = Endmember3(x, cx, y, cy, z, cz)

function Endmember(;x::Datum=Constant(NaN), cx::Datum=Constant(NaN), y::Datum=Constant(NaN), cy::Datum=Constant(NaN), z::Datum=Constant(NaN), cz::Datum=Constant(NaN)) 
    
    @assert Constant(NaN) ∉ (x,cx) "Insufficient components declared."
    @assert iseven(count(x -> x==Constant(NaN),(x,y,z,cx,cy,cz))) "All components must have paired isotope composition and concentration data."

    out = if z == cz == y == cy == Constant(NaN)
        Endmember1(x,cx)
    elseif z == cz == Constant(NaN)
        Endmember2(x, cx, y, cy)
    else
        Endmember3(x, cx, y, cy, z, cz)
    end
    out
end
"""

    Endmember1 <: Endmember 

1-dimensional `Endmember` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x

see also: [`Endmember`](@ref), [`Datum`](@ref)

"""
struct Endmember1 <: Endmember
    x::Datum
    cx::Datum
end

"""

    Endmember2 <: Endmember 

2-dimensional `Endmember` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y

see also: [`Endmember`](@ref), [`Datum`](@ref)

"""
struct Endmember2 <: Endmember
    x::Datum
    cx::Datum
    y::Datum
    cy::Datum
end

"""

    Endmember3 <: Endmember 

3-dimensional `Endmember` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y
`z` | Isotopic composition of z
`cz`| Concentration of z

see also: [`Endmember`](@ref), [`Datum`](@ref)

"""
struct Endmember3 <: Endmember
    x::Datum
    cx::Datum
    y::Datum
    cy::Datum
    z::Datum
    cz::Datum
end

## Prior

"""
    Prior <: Data

Abstract supertype for `Prior_` instances, where the `_` indicates dimensionality. Represents a suite of [`Endmember`] data.

see also: [`Endmember`](@ref), [`Prior2`](@ref), [`Prior3`](@ref)

---

    Prior(A::Endmember, B::Endmember)

Returns a `Prior2`

    Prior(A::Endmember, B::Endmember, C::Endmember)

Returns a `Prior3`

Note: each `Endmember` must be of the same subtype (i.e. consistent number of components within each endmember of a system).

"""
abstract type Prior{T<:Endmember} <: Data end
Prior(A::T, B::T) where T<:Endmember = Prior2(A,B)
Prior(A::T, B::T, C::T) where T<:Endmember = Prior3(A,B,C)

"""
    Prior2 <: Prior

2-dimensional `Prior` instance.

|Fields ||
|:--|:--|
`A` | Endmember/component A
`B` | Endmember/component B

see also: [`Prior`](@ref)
"""
struct Prior2{T<:Endmember} <: Prior{T}
    A::T
    B::T
end

"""
    Prior3 <: Prior

3-dimensional `Prior` instance.

|Fields ||
|:--|:--|
`A` | Endmember/component A
`B` | Endmember/component B
`C` | Endmember/component C

see also: [`Prior`](@ref)

"""
struct Prior3{T<:Endmember} <: Prior{T}
    A::T
    B::T
    C::T
end

## Fraction

"""
    Fraction <: IsoMixType

Abstract supertype for `Fraction_` instances that contain vectors of fractional combinations of `_` components/endmembers.

see also: [`Fraction2`](@ref), [`Fraction3`](@ref)

---

    Fraction(Amin, Amax; n=101)

Returns a `Fraction2` with `n` linearly spaced fractional combinations of components `A` and `B` that each sum to 1, given minimum and maximum A fractions of `Amin` and `Amax`, respectively. Since`n` gives the linear fractions of each component, this constructor returns `n`-length vectors of each component. Both `Amin` and `Amax` must fall within [0,1].

# Example

    julia> Fraction(0,1,n=3)
    Fraction2([0.0, 0.5, 1.0], [1.0, 0.5, 0.0], 3)
---

    Fraction(Amin, Amax, Bmin, Bmax; n=101)

Returns a `Fraction3` with linearly spaced fractional combinations of components `A`, `B`, and `C` that each sum to 1, given respective minima and maxima of `A` (`Amin` and `Amax`) and `B` (`Bmin` and `Bmax`). Since `n` gives the linear fractions of each component, this constructor returns a `n(n÷2+1)`-length vector for each component. `Amin`,`Amax`, `Bmin`, and `Bmax` must fall within [0,1].

# Example

    julia> f=Fraction(0,1,0,1,n=3)
    Fraction3([0.0, 0.0, 0.0, 0.5, 0.5, 1.0], [0.0, 0.5, 1.0, 0.0, 0.5, 0.0], [1.0, 0.5, 0.0, 0.5, 0.0, 0.0], 6)

---

    Fraction(dim, n=101)

Given a provided number of components/dimensions `dim` (2 or 3), returns a corresponding `Fraction` instance with component concentrations of 0 to 1. 

# Example
    julia> Fraction(2,n=3)
    Fraction2([0.0, 0.5, 1.0], [1.0, 0.5, 0.0], 3)

"""
abstract type Fraction <: IsoMixType end 
function Fraction(D::Int=2; n::Int=101)
    if D==2
        Fraction(0,1,n=n)
    elseif D==3
        Fraction(0,1,0,1, n=n)
    else
        error("Must D must be a value of either 2 or 3.")
    end
end

function Fraction(Amin::Number, Amax::Number; n::Int=101)
    @assert 0 ≤ Amin < Amax ≤ 1
    
    fA = collect(range(float(Amin), float(Amax), n))
    fB = 1.0 .- fA
    Fraction2(fA,fB,n)
end

function Fraction(Amin::Number, Amax::Number, Bmin::Number, Bmax::Number; n::Int=101)

    @assert 0 ≤ Amin < Amax ≤ 1
    @assert 0 ≤ Bmin < Bmax ≤ 1

    Arange = range(float(Amin), float(Amax), n)
    Brange = range(float(Bmin), float(Bmax), n)

    nv = n * (div(n,2)+1)
    fA = Vector{eltype(Arange)}(undef,nv)
    fB, fC = similar(fA), similar(fA)

    vi = 0 
    for i in eachindex(Arange), j in eachindex(Brange)
        fAi, fBi =  Arange[i], Brange[j]
        fAB = fAi + fBi
        fCi = 1.0 - fAB 
        sumto1 = fAB <= 1.0

        vi += ifelse(sumto1, 1, 0)
        fA[vi] = ifelse(sumto1, fAi, fA[vi])
        fB[vi] = ifelse(sumto1, fBi, fB[vi])
        fC[vi] = ifelse(sumto1, fCi, fC[vi])
    end
    Fraction3(fA,fB,fC,nv)
end


"""
    Fraction2 <: Fraction

2-dimensional `Fraction` instance.

|Fields ||
|:--|:--|
`A` | Endmember/component A
`B` | Endmember/component B

see also: [`Fraction`](@ref)

"""
struct Fraction2 <: Fraction
    A::Vector{Float64}
    B::Vector{Float64}
    n::Int
end

"""
    Fraction3 <: Fraction

3-dimensional `Fraction` instance.

|Fields ||
|:--|:--|
`A` | Endmember/component A
`B` | Endmember/component B
`C` | Endmember/component C

see also: [`Fraction`](@ref)

"""
struct Fraction3 <: Fraction
    A::Vector{Float64}
    B::Vector{Float64}
    C::Vector{Float64}
    n::Int
end

## Model

"""

    Model <: IsoMixType

Abstract supertype for `Model_` instances, where the `_` indicates dimensionality.

see also: [`Model1`](@ref), [`Model2`](@ref), [`Model3`](@ref)

---

    Model(f::Fraction, s::SysDraw)

Generate a `Model` instance 

"""
abstract type Model <: IsoMixType end

Model(f::Fraction,::SysDraw{EmDraw1})= Model1(Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n))

Model(f::Fraction,::SysDraw{EmDraw2}) = Model2(Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n),Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n))

Model(f::Fraction,::SysDraw{EmDraw3}) = Model3(Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n),Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n),Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n))

"""

    Model1 <: Model

1-dimensional `Model` instance: a single species tracking isotopic composition and concentration.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x

see also: [`Model`](@ref)

"""
struct Model1 <: Model
    x::Vector{Float64}
    cx::Vector{Float64}
end

"""

    Model2 <: Model

2-dimensional `Model` instance: two species tracking isotopic composition and concentration.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y

see also: [`Model`](@ref)

"""
struct Model2 <: Model
    x::Vector{Float64}
    cx::Vector{Float64}
    y::Vector{Float64}
    cy::Vector{Float64}
end

"""

    Model3 <: Model

3-dimensional `Model` instance: three species tracking isotopic composition and concentration.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y
`z` | Isotopic composition of z
`cz`| Concentration of z

see also: [`Model`](@ref)

"""
struct Model3 <: Model
    x::Vector{Float64}
    cx::Vector{Float64}
    y::Vector{Float64}
    cy::Vector{Float64}
    z::Vector{Float64}
    cz::Vector{Float64}
end

## Measurements

"""

    Measurements <: Data

Abstract supertype for Measurements_ instances, where _ indicates dimensionality. Represents a dataset of measured (normally distributed) data as vectors of [`Norm`](@ref) instances. Absent `Measurement` fields are represented as vectors of [`Unconstrained`](@ref) and do not affect [`loglikelihood`](@ref) calculations. 

see also: [`Measurements`](@ref), [`Measurements1`](@ref), [`Measurements2`](@ref), [`Measurements3`](@ref)

# Construction

The constructor function `Measurements(...)` takes `nx2` matrices, where `n` represents the number of measured samples in the dataset, and `n` must be the same for all provided matrices. Each matrix represents a series of measurements (e.g. isotopic composition or abundance measurement) for a given species. The table below defines all possible fields of a `Measurements` and which `Measurements` subtypes contain them (denoted with a ✓):

|Fields || `Measurements1` | `Measurements2` | `Measurements3`
|:--|:--| :--: | :-: | :-: |
`x` | Isotopic composition of x | ✓ | ✓ | ✓ |
`cx`| Concentration of x | ✓ | ✓ | ✓ |
|||||
`y` | Isotopic composition of y |  | ✓ | ✓ |
`cy`| Concentration of y |  | ✓ | ✓ |
|||||
`z` | Isotopic composition of z |  |  | ✓ |
`cz`| Concentration of z |  |  | ✓ |

Note that each row among the `Measurement` fields corresponds to a specific sample, so provided matrices must consistently reflect the same sample in each row. Measurements with missing values must be replaced with `NaN`. 

## Inputing data into `Measurements()`

The constructor function accepts a list of `Measurements` instances:

Without concentration measurements:

    Measurements(x, y) -> Measurements2 
    Measurements(x, y, z) -> Measurements3
    
With concentration measurements:

    Measurements(x=..., cx=...) -> Measurements1
    Measurements(x, cx, y, cy) -> Measurements2
    Measurements(x, cx, y, cy, z, cz) -> Measurements3

For all other datasets configurations, declare measurements with keywords:

    Measurements(; x, cx, y, cy, z, cz)

...which determines the appropriate subtype and all undeclared fields return vectors of `Unconstrained` instances.

## Examples

```julia

julia> Measurements(x = [1 2; 3 4], cx = [5 6; 7 8])
Measurements1(Datum[Norm(1.0, 2.0), Norm(3.0, 4.0)], Datum[Norm(5.0, 6.0), Norm(7.0, 8.0)])

julia> m2d = Measurements([1 2; 3 4], [5 6; 7 NaN]); 

julia> m2d.x
2-element Vector{Datum}:
 Norm(1.0, 2.0)
 Norm(3.0, 4.0)

julia> m2d.cx
2-element Vector{Datum}:
 Unconstrained()
 Unconstrained()

julia> m2d.y
2-element Vector{Datum}:
 Norm(5.0, 6.0)
 Norm(7.0, NaN)

julia> m2d = Measurements([1 2; 3 4], [5 6; 7 8; 9 10])
ERROR: AssertionError: x & y must have the same number of entries (rows). Use NaN for missing values.

```

"""
abstract type Measurements <: Data end

function Measurements(x::Matrix, y::Matrix)
    @assert size(x,1) == size(y,1) "x & y must have the same number of entries (rows). Use NaN for missing values."
    @assert size(x,2) == size(y,2) == 2 "all entries must have only 2 columns"
    
    n = size(x,1)
    xd, yd = (Vector{Norm}(undef,n) for i=1:2)
    
    @inbounds for i in 1:n
        xd[i] = Norm(x[i,1], x[i,2])
        yd[i] = Norm(y[i,1], y[i,2])
    end
    c = fill(Unconstrained(),n)
    Measurements2(xd,c,yd,c)
end

function Measurements(x::Matrix, y::Matrix, z::Matrix)
    @assert size(x,1) == size(y,1) == size(z,1) "x, y, & z must have the same number of entries (rows). Use NaN for missing values."
    @assert size(x,2) == size(y,2) == size(z,2) == 2 "all entries must have only 2 columns"
    
    n = size(x,1)
    xd,yd,zd = (Vector{Norm}(undef,n) for i=1:3)
    
    @inbounds for i in 1:n
        xd[i] = Norm(x[i,1], x[i,2])
        yd[i] = Norm(y[i,1], y[i,2])
        zd[i] = Norm(z[i,1], z[i,2])
    end
    c = fill(Unconstrained(),n)
    Measurements3(xd,c,yd,c,zd,c)
end

function Measurements(x::Matrix, cx::Matrix, y::Matrix, cy::Matrix)
    @assert size(x,1) == size(cx,1) == size(y,1) == size(cy,1) "x, cx, y, & cy must have the same number of entries (rows). Use NaN for missing values."
    @assert size(x,2) == size(cx,2) == size(y,2) == size(cy,2) == 2 "all entries must have only 2 columns"
    
    n = size(x,1)
    xd,cxd,yd,cyd = (Vector{Norm}(undef,n) for i=1:4)
    
    @inbounds for i in 1:n
        xd[i] = Norm(x[i,1], x[i,2])
        yd[i] = Norm(y[i,1], y[i,2])
        cxd[i] = Norm(cx[i,1], cx[i,2])
        cyd[i] = Norm(cy[i,1], cy[i,2])
    end
    Measurements2(xd,cxd,yd,cyd)
end

function Measurements(x::Matrix, cx::Matrix, y::Matrix, cy::Matrix, z::Matrix, cz::Matrix)
    @assert size(x,1) == size(cx,1) == size(y,1) == size(cy,1) == size(z,1) == size(cz,1) "x, cx, y, cy, z, & cz must have the same number of entries (rows). Use NaN for missing values."
    @assert size(x,2) == size(cx,2) == size(y,2) == size(cy,2) == size(z,2) == size(cz,2) == 2 "all entries must have only 2 columns"

    n = size(x,1)
    xd,cxd,yd,cyd,zd,czd = (Vector{Norm}(undef,n) for i=1:6)

    @inbounds for i in 1:n
        xd[i] = Norm(x[i,1], x[i,2])
        yd[i] = Norm(y[i,1], y[i,2])
        zd[i] = Norm(z[i,1], z[i,2])
        cxd[i] = Norm(cx[i,1], cx[i,2])
        cyd[i] = Norm(cy[i,1], cy[i,2])
        czd[i] = Norm(cz[i,1], cz[i,2])
    end
    c = fill(Unconstrained(),n)
    Measurements3(xd,c,yd,c,zd,czd)
end

function Measurements(; x::Matrix=[Inf Inf], cx::Matrix=[Inf Inf], y::Matrix=[Inf Inf], cy::Matrix=[Inf Inf], z::Matrix=[Inf Inf], cz::Matrix=[Inf Inf])
    @assert size(x,2) == size(cx,2) == size(y,2) == size(cy,2) == size(z,2) == size(cz,2) == 2 "All entries must have only 2 columns"
    @assert x != [Inf Inf] "x must be defined. Provide a matrix of measurements for x."
    
    n = size(x,1)
    null = [Inf Inf] 
    sizetest = cx == null || size(cx,1) == n 
    sizetest *= y == null  || size(y,1) == n 
    sizetest *= cy == null || size(cy,1) == n
    sizetest *= z == null  || size(z,1) == n
    sizetest *= cz == null || size(cz,1) == n 
    @assert sizetest "All provided inputs must have the same number of entries (rows). Use NaN for missing values."
    
    xd = Vector{Norm}(undef,n)
    cxd = cx == null ? Vector{Unconstrained}(undef,n) : Vector{Norm}(undef,n)
    yd = y == null ? Vector{Unconstrained}(undef,n) : Vector{Norm}(undef,n)
    cyd = cy == null ? Vector{Unconstrained}(undef,n) : Vector{Norm}(undef,n)
    zd = z == null ? Vector{Unconstrained}(undef,n) : Vector{Norm}(undef,n)
    czd = cz == null ? Vector{Unconstrained}(undef,n) : Vector{Norm}(undef,n)
    
    @inbounds for i in 1:n
        xd[i] = Norm(x[i,1], x[i,2])
        cx == null || (cxd[i] = Norm(cx[i,1], cx[i,2]))
        y == null || (yd[i] = Norm(y[i,1], y[i,2]))
        cy == null || (cyd[i] = Norm(cy[i,1], cy[i,2]))
        z == null || (zd[i] = Norm(z[i,1], z[i,2]))
        cz == null || (czd[i] = Norm(cz[i,1], cz[i,2]))
    end
    
    if eltype(czd) == Norm || eltype(zd) == Norm
        Measurements3(xd,cxd,yd,cyd,zd,czd)
    elseif eltype(cyd) == Norm || eltype(yd) == Norm
        Measurements2(xd,cxd,yd,cyd)
    else
        Measurements1(xd,cxd)
    end
end


"""
    
    Measurements1 <: Measurements

[`Measurements`](@ref) instance reflecting isotope and abundance measurements of 1 species: x.

"""
struct Measurements1 <: Measurements
    x::Vector{Datum} 
    cx::Vector{Datum} 
end

"""

    Measurements2 <: Measurements

[`Measurements`](@ref) instance reflecting isotope and abundance measurements of 2 species: x, y.

"""
struct Measurements2 <: Measurements
    x::Vector{Datum} 
    cx::Vector{Datum} 
    y::Vector{Datum} 
    cy::Vector{Datum} 
end

"""

    Measurements3 <: Measurements

[`Measurements`](@ref) instance reflecting isotope and abundance measurements of 3 species: x, y, z.

"""
struct Measurements3 <: Measurements
    x::Vector{Datum} 
    cx::Vector{Datum} 
    y::Vector{Datum} 
    cy::Vector{Datum} 
    z::Vector{Datum} 
    cz::Vector{Datum} 
end