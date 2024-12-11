
## Fraction

"""
    Fraction <: Calculations

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
abstract type Fraction <: Calculations end 
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
    @assert n >= 3
    fA = collect(range(float(Amin), float(Amax), n))
    fB = 1.0 .- fA
    Fraction2(fA,fB,n)
end

function Fraction(Amin::Number, Amax::Number, Bmin::Number, Bmax::Number; n::Int=101)

    @assert 0 ≤ Amin < Amax ≤ 1
    @assert 0 ≤ Bmin < Bmax ≤ 1
    @assert n >= 3
    
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

        vi += sumto1 
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

    Model <: Calculations

Abstract supertype for `Model_` instances, where the `_` indicates dimensionality.

see also: [`Model1`](@ref), [`Model2`](@ref), [`Model3`](@ref)

---

    Model(f::Fraction, s::SysDraw)

Generate a `Model` instance 

"""
abstract type Model <: Calculations end

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

## Chains

"""

    EmChains <: Calculations

Short for `EndmemberChains`. Abstract supertype for `EmC_` instances, where the `_` indicates dimensionality. Represents a Markov chain posterior for the constituent compositions of an endmember within a natural system. 

see also: [`Chains`](@ref), [`EmC1`](@ref), [`EmC2`](@ref), [`EmC3`](@ref)

### Construction

    EmChains(n, <:Endmember)
    
To construct, seed with a number of chain iterations `n` and an `Endmember` instance to seed dimensionality, e.g.

`EmChains(n,::Endmember2)` -> `EmC2`

"""
abstract type EmChains <: Calculations end

EmChains(n::Int,::Endmember1) = EmC1(Vector{Float64}(undef,n),Vector{Float64}(undef,n))

EmChains(n::Int,::Endmember2) = EmC2(Vector{Float64}(undef,n), Vector{Float64}(undef,n), Vector{Float64}(undef,n), Vector{Float64}(undef,n))

EmChains(n::Int,::Endmember3) =  EmC3(Vector{Float64}(undef,n), Vector{Float64}(undef,n), Vector{Float64}(undef,n), Vector{Float64}(undef,n), Vector{Float64}(undef,n), Vector{Float64}(undef,n))

"""

    EmC1 <: EmChains 

1-dimensional `EmChains` instance.

|Fields ||
|:--|:--|
`x` | Markov chain of isotopic compositions of x
`cx`| Markov chain of concentrations of x

see also: [`EmChains`](@ref), [`Chains`](@ref), [`Endmember`](@ref)

"""
struct EmC1 <: EmChains
    x::Vector{Float64}
    cx::Vector{Float64}
end

"""

    EmC2 <: EmChains 

2-dimensional `EmChains` instance.

|Fields ||
|:--|:--|
`x` | Markov chain of isotopic compositions of x
`cx`| Markov chain of concentrations of x
`y` | Markov chain of isotopic compositions of y
`cy`| Markov chain of concentrations of y

see also: [`EmChains`](@ref), [`Chains`](@ref), [`Endmember`](@ref)

"""
struct EmC2 <: EmChains
    x::Vector{Float64}
    cx::Vector{Float64}
    y::Vector{Float64}
    cy::Vector{Float64}
end

"""

    EmC3 <: EmChains 

3-dimensional `EmChains` instance.

|Fields ||
|:--|:--|
`x` | Markov chain of isotopic compositions of x
`cx`| Markov chain of concentrations of x
`y` | Markov chain of isotopic compositions of y
`cy`| Markov chain of concentrations of y
`z` | Markov chain of isotopic compositions of z
`cz`| Markov chain of concentrations of z

see also: [`EmChains`](@ref), [`Chains`](@ref), [`Endmember`](@ref)

"""
struct EmC3 <: EmChains
    x::Vector{Float64}
    cx::Vector{Float64}
    y::Vector{Float64}
    cy::Vector{Float64}
    z::Vector{Float64}
    cz::Vector{Float64}
end



"""





---

    Prior(A<:EmChains, B::Endmember)

Returns a `Prior2`

    Prior(A::Endmember, B::Endmember, C::Endmember)

Returns a `Prior3`

Note: each `Endmember` must be of the same subtype (i.e. consistent number of components within each endmember of a system).

"""

"""

    Chains <: Calculations

Abstract supertype for `Chains_` instances, where the `_` indicates dimensionality. Represents a suite of EndmemberChains ([`EmChains`](@ref)). 

see also: [`EmChains`](@ref), [`Chains2`](@ref), [`Chains3`](@ref), [`Prior`](@ref)

### Construction

    Chains(n, <:Prior)
    
To construct, seed with a number of chain iterations `n` and a Prior instance to seed dimensionality, e.g.

`Chains(n,::Prior2{Endmember3})` -> `Chains2{EmC3}`

"""
abstract type Chains{T<:EmChains} <: Calculations end

Chains(n::Int,p::Prior2) =  Chain2( EmChains(n, p.A), EmChains(n, p.A), Vector{Float64}(undef,n), BitVector(undef,n) )

Chains(n::Int,p::Prior3) =  Chain3( EmChains(n, p.A), EmChains(n, p.A), EmChains(n, p.A), Vector{Float64}(undef,n), BitVector(undef,n) )

"""
    Chains2 <: Prior

2-dimensional `Chains` instance.

|Fields ||
|:--|:--|
`A` | Markov chains of endmember/component A
`B` | Markov chains of endmember/component B

see also: [`Chains`](@ref)

"""
struct Chains2{T<:EmChains} <: Chains
    A::T
    B::T
    ll::Vector{Float64}
    accept::BitVector
end

"""
    Chains3 <: Prior

3-dimensional `Chains` instance.

|Fields ||
|:--|:--|
`A` | Markov chains of endmember/component A
`B` | Markov chains of endmember/component B
`C` | Markov chains of endmember/component C

see also: [`Chains`](@ref)

"""
struct Chains3{T<:EmChains} <: Chains
    A::T
    B::T
    C::T
    ll::Vector{Float64}
    accept::BitVector
end

