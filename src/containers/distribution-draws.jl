
## EmDraw 

"""
    EmDraw <: DistributionDraws

Short for `EndmemberDraw`. Abstract supertype for `EmDraw_` instances, where the `_` indicates dimensionality. Represents the composition of an endmember within a natural system. 

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

## SysDraw

"""

    SysDraw <: Data

Short for `SystemDraw`. Abstract supertype for `SysDraw_` instances, where the `_` indicates dimensionality. Represents a model system composition charaterized by its `EmDraw` fields.

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