"""
    IsoMixType

Supertype containing all custom types and structs in the IsoMix package:

IsoMixType

├─ [`Component`](@ref)

│  ├─ Component1

│  ├─ Component2

│  └─ Component3

├─ [`Data`](@ref)

│  ├─ Data1

│  ├─ Data2

│  └─ Data3

├─ [`Datum`](@ref)

│  ├─ [`Constant`](@ref)

│  ├─ [`Norm`](@ref)

│  ├─ [`Unf`](@ref)

│  └─ [`logNorm`](@ref)

├─ [`Fraction`](@ref)

│  ├─ Fraction2

│  └─ Fraction3

├─ [`Model`](@ref)

│  ├─ Model1

│  ├─ Model2

│  └─ Model3

├─ [`Prior`](@ref)

│  ├─ Prior2

│  └─ Prior3

└─ [`System`](@ref)

├─ System2

└─ System3


# Example

`julia> using AbstractTrees,IsoMix; AbstractTrees.children(d) = subtypes(d); print_tree(IsoMixType)`

"""
abstract type IsoMixType end

## Component

"""
    Component <: IsoMixType

Abstract supertype for `Component_` instances, where the `_` indicates dimensionality. Represents the composition of a component or endmember within a natural system. 

see also: [`Component1`](@ref), [`Component2`](@ref), [`Component3`](@ref)

---

`Component(x, cx)` -> `Component1`

`Component(x, cx, y, cy)` -> `Component2`

`Component(x, cx, y, cy, z, cz)` -> `Component3`

The constructor function accepts values as a list of arguments or as keyword assignments, e.g.: `Component(x=1, cx=2)`

"""
abstract type Component <: IsoMixType end
Component(x::Number, cx::Number) = Component1(x, cx)

Component(x::Number, cx::Number, y::Number, cy::Number) = Component2(x, cx, y, cy)

Component(x::Number, cx::Number, y::Number, cy::Number, z::Number, cz::Number) = Component3(x, cx, y, cy, z, cz)

function Component(;x::Number=NaN, cx::Number=NaN, y::Number=NaN, cy::Number=NaN, z::Number=NaN, cz::Number=NaN) 
    out = if (z == z) & (cz == cz)
        Component3(x, cx, y, cy, z, cz)
    elseif (y==y) & (cy==cy)
        Component2(x, cx, y, cy)
    elseif (x==x) & (cx==cx)
        Component1(x,cx)
    else
        error("Improper Component definition")
    end
    out
end

"""
    Component1 <: Component

1-dimensional `Component` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x

see also: [`Component`](@ref)

"""
struct Component1 <: Component
    x::Float64
    cx::Float64
end 

"""
    Component2 <: Component

2-dimensional `Component` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y

see also: [`Component`](@ref)

"""
struct Component2 <: Component
    x::Float64
    cx::Float64
    y::Float64
    cy::Float64
end 

"""
    Component3 <: Component

3-dimensional `Component` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y
`z` | Isotopic composition of z
`cz`| Concentration of z

see also: [`Component`](@ref)

"""
struct Component3 <: Component
    x::Float64
    cx::Float64
    y::Float64
    cy::Float64
    z::Float64
    cz::Float64
end



"""

    IsoMix.countcomponents(c::Component)

Count the number of elemental dimensions (1-3) of `c`.

# Example

    julia> IsoMix.countcomponents(Component(1,2,3,4,5,6))
    3

"""
countcomponents(::C) where C <: Component = div(fieldcount(C),2)

## System

"""

    Prior <: IsoMixType

Abstract supertype for `System_` instances, where the `_` indicates dimensionality. Represents a model system composition charaterized its `Component` fields.

see also: [`Prior2`](@ref), [`Prior3`](@ref)

---
The constructor function accepts values as a list of arguments:

    System(A::Component, B::Component)
Returns a `System2`

    System(A::Component, B::Component, C::Component)
Returns a `System3`

... or as keyword assignments, e.g.: 
    System(A=..., B=...)

"""
abstract type System{T} <: IsoMixType end
System(A::T, B::T) where T <: Component= System2(A,B)
System(A::T, B::T, C::T) where T <: Component = System3(A,B,C)

"""
    System2 <: System

2-dimensional `System` instance.

|Fields ||
|:--|:--|
`A` | Endmember/component A
`B` | Endmember/component B

see also: [`System`](@ref)
"""
@kwdef struct System2{T<:Component} <: System{T}
    A::T
    B::T
end

"""
    System3 <: System

3-dimensional `System` instance.

|Fields ||
|:--|:--|
`A` | Endmember/component A
`B` | Endmember/component B
`C` | Endmember/component C

see also: [`System`](@ref)
"""
@kwdef struct System3{T<:Component} <: System{T}
    A::T
    B::T
    C::T
end

## Datum

"""
    Datum <: IsoMixType

Abstract supertype for `Datum` instances.

see also: [`Constant`](@ref), [`Unf`](@ref), [`Norm`](@ref), [`logNorm`](@ref)

"""
abstract type Datum <: IsoMixType end 

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

## Data

"""

    Data <: IsoMixType

Abstract supertype for `Data_` instances, where the `_` indicates dimensionality.

see also: [`Data1`](@ref), [`Data2`](@ref), [`Data3`](@ref)

"""
abstract type Data <: IsoMixType end

Data(x::Datum, cx::Datum) = Data1(x, cx)

Data(x::Datum, cx::Datum, y::Datum, cy::Datum) = Data2(x, cx, y, cy)

Data(x::Datum, cx::Datum, y::Datum, cy::Datum, z::Datum, cz::Datum) = Data3(x, cx, y, cy, z, cz)

function Data(;x::Datum=Constant(NaN), cx::Datum=Constant(NaN), y::Datum=Constant(NaN), cy::Datum=Constant(NaN), z::Datum=Constant(NaN), cz::Datum=Constant(NaN)) 
    
    @assert Constant(NaN) ∉ (x,cx) "Insufficient components declared."
    @assert iseven(count(x -> x==Constant(NaN),(x,y,z,cx,cy,cz))) "All components must have paired isotope composition and concentration data."

    out = if z == cz == y == cy == Constant(NaN)
        Data1(x,cx)
    elseif z == cz == Constant(NaN)
        Data2(x, cx, y, cy)
    else
        Data3(x, cx, y, cy, z, cz)
    end
    out
end
"""

    Data1 <: Data 

1-dimensional `Data` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x

see also: [`Datum`](@ref)

"""
struct Data1 <: Data
    x::Datum
    cx::Datum
end

"""

    Data2 <: Data 

2-dimensional `Data` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y

see also: [`Datum`](@ref)

"""
struct Data2 <: Data
    x::Datum
    cx::Datum
    y::Datum
    cy::Datum
end

"""

    Data3 <: Data 

3-dimensional `Data` instance.

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y
`z` | Isotopic composition of z
`cz`| Concentration of z

see also: [`Datum`](@ref)

"""
struct Data3 <: Data
    x::Datum
    cx::Datum
    y::Datum
    cy::Datum
    z::Datum
    cz::Datum
end

## Prior

"""
    Prior <: IsoMixType

Abstract supertype for `Prior_` instances, where the `_` indicates dimensionality. Represents a suite of data that provide a prior estimate for the endmember compositions of a natural system.

see also: [`Prior2`](@ref), [`Prior3`](@ref)

---

    Prior(A::Data, B::Data)

Returns a `Prior2`

    Prior(A::Data, B::Data, C::Data)

Returns a `Prior3`

Note: each `Data` must be of the same subtype (i.e. consistent number of components within each endmember of a system).

"""
abstract type Prior <: IsoMixType end
Prior(A::Data, B::Data) = Prior2(A,B)
Prior(A::Data, B::Data, C::Data) = Prior3(A,B,C)

"""
    Prior2 <: Prior

2-dimensional `Prior` instance.

|Fields ||
|:--|:--|
`A` | Endmember/component A
`B` | Endmember/component B

see also: [`Prior`](@ref)
"""
struct Prior2{T} <:Prior where T <: Data
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
struct Prior3{T} <:Prior where T<:Data
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

    Fraction(D, n=101)

Given a provided number of components `D` (2 or 3), returns a corresponding `Fraction` instance with component concentrations of 0 to 1. 

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

    Model(f::Fraction, s::System)

Generate a `Model` instance 

"""
abstract type Model <: IsoMixType end

Model(f::Fraction,::System{Component1})= Model1(Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n))

Model(f::Fraction,::System{Component2}) = Model2(Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n),Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n))

Model(f::Fraction,::System{Component3}) = Model3(Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n),Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n),Vector{Float64}(undef,f.n), Vector{Float64}(undef,f.n))

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

"""

    Measurements <: IsoMixType

Struct containing two vector fields corresponding to the mean `m` and 1σ standard deviation (or standard error) `s` of measurements with normally distributed uncertainties (typical of geochemical analyses). All `IsoMix` applications of this struct require `length(m)==length(s)`. Measurements with missing values may either be removed or replaced with `NaN`. 

see also: [`DataSet`](@ref)

---

    Measurements(m,s) 
    Measurements(; m, s)
    Measurements(m; s)
    Measurements() --> returns a null instance of Measurements([],[])

# Examples 

    julia> Measurements([2.,3.],[0.5,0.5])
    Measurements([2.0, 3.0], [0.5, 0.5])

    julia> Measurements(m=[2.,3.],s=[0.5,0.5])
    Measurements([2.0, 3.0], [0.5, 0.5])    

    julia> Measurements([2.,3.],s=[0.5,0.5])
    Measurements([2.0, 3.0], [0.5, 0.5])

    julia> Measurements()
    Measurements(Float64[], Float64[])


"""
struct Measurements <: IsoMixType
    m::Vector{Float64}
    s::Vector{Float64}
end
Measurements(; m::Vector{<:Number}=Float64[], s::Vector{<:Number}=Float64[]) = Measurements(m,s)
Measurements(m::Vector{<:Number}; s::Vector{<:Number}=Float64[]) = Measurements(m,s)



"""

    DataSet <: IsoMixType

  Abstract supertype for DataSet_ instances, where _ indicates dimensionality. Represents a dataset of measured (normally distributed) data as `Measurements`.

  see also: [`Measurements`](@ref), [`DataSet1`](@ref), [`DataSet2`](@ref), [`DataSet3`](@ref)

  ---

The constructor function accepts a list of `Measurements` instances:

  DataSet(x, cy) -> DataSet1
  DataSet(x, cx, y, cy) -> DataSet2
  DataSet(x, cx, y, cy, z, cz) -> DataSet3

For inconsistent datasets, you may declare extant data with keywords:

    DataSet(; x, cx, y, cy, z, cz)

`DataSet` determines the appropriate subtype and all undeclared fields are empty: `Measurements([],[])`

NOTE: the constructor requires that `length(q.m)==length(q.s)` for each `Measurements` instance `q`.

"""
abstract type DataSet <: IsoMixType end
function DataSet(x::T, cx::T) where T<:Measurements
    @assert length(x.m) == length(x.s)
    @assert length(cx.m) == length(cx.s)
    DataSet1(x,cx)
end
function DataSet(x::T, cx::T, y::T, cy::T) where T<:Measurements
    @assert length(x.m) == length(x.s)
    @assert length(cx.m) == length(cx.s)
    @assert length(y.m) == length(y.s)
    @assert length(cy.m) == length(cy.s)
    DataSet2(x,cx,y,cy)
end
function DataSet(x::T, cx::T, y::T, cy::T, z::T, cz::T) where T<:Measurements
    @assert length(x.m) == length(x.s)
    @assert length(cx.m) == length(cx.s)
    @assert length(y.m) == length(y.s)
    @assert length(cy.m) == length(cy.s)
    @assert length(z.m) == length(z.s)
    @assert length(cz.m) == length(cz.s)
    DataSet3(x,cx,y,cy,z,cz)
end
function DataSet(; x::T=Measurements(), cx::T=Measurements(), y::T=Measurements(), cy::T=Measurements(), z::T=Measurements(), cz::T=Measurements()) where T<: Measurements
    if (!isempty(z.m) | !isempty(cz.m))
        DataSet(x,cx,y,cy,z,cz)
    elseif (!isempty(y.m) | !isempty(cy.m))
        DataSet(x,cx,y,cy)
    else
        DataSet(x,cx)
    end
end

"""

    DataSet1 <: DataSet

`DataSet` instance reflecting measurements of 1 species:

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x

see also: [`DataSet`](@ref)

"""
struct DataSet1 <: DataSet
    x::Measurements
    cx::Measurements
end

"""

    DataSet2 <: DataSet

`DataSet` instance reflecting measurements of 2 species:

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y

see also: [`DataSet`](@ref)

"""
struct DataSet2 <: DataSet
    x::Measurements
    cx::Measurements
    y::Measurements
    cy::Measurements
end

"""

    DataSet3 <: DataSet

`DataSet` instance reflecting measurements of 3 species:

|Fields ||
|:--|:--|
`x` | Isotopic composition of x
`cx`| Concentration of x
`y` | Isotopic composition of y
`cy`| Concentration of y
`z` | Isotopic composition of z
`cz`| Concentration of z

see also: [`DataSet`](@ref)

"""
struct DataSet3 <: DataSet
    x::Measurements
    cx::Measurements
    y::Measurements
    cy::Measurements
    z::Measurements
    cz::Measurements
end