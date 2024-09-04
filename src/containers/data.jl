

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