

struct Component
    x::Float64
    Cx::Float64
    y::Float64
    Cy::Float64
    z::Float64
    Cz::Float64
end

Component(x::Number, Cx::Number) = Component(x, Cx, NaN, NaN, NaN, NaN)
Component(x::Number, Cx::Number, y::Number, Cy::Number) = Component(x, Cx, y, Cy, NaN, NaN)

function countcomponents(c::Component)
    f = fieldnames(Component)
    n = 0
    @inbounds for i in f
        n += ifelse(isnan(getfield(c,i)),0,1)
    end
    iseven(n) || @warn "Non-even non-NaN components. Inspect input."
    n
end


abstract type System end

struct System2 <: System
    A::Component
    B::Component
end

struct System3 <: System
    A::Component
    B::Component
    C::Component
end

System(A::T, B::T) where T <: Component= System2(A,B)
System(A::T, B::T, C::T) where T <: Component = System3(A,B,C)


abstract type Datum end 

struct Norm <: Datum
    mu::Float64
    sig::Float64
end

struct logNorm <: Datum
    mu::Float64
    sig::Float64
end

struct Unf <: Datum
    a::Float64
    b::Float64
end

struct Constant <: Datum
    x::Float64
end

abstract type Data end

struct Data2 <: Data
    x::Datum
    Cx::Datum
    y::Datum
    Cy::Datum
end

struct Data3 <: Data
    x::Datum
    Cx::Datum
    y::Datum
    Cy::Datum
    z::Datum
    Cz::Datum
end

abstract type Prior end

struct Prior2{T} <:Prior where T <: Data
    A::T
    B::T
end

struct Prior3{T} <:Prior where T <: Data
    A::T
    B::T
    C::T
end

Prior(A::Data, B::Data) = Prior2(A,B)
Prior(A::Data, B::Data, C::Data) = Prior3(A,B,C)

abstract type Fraction end 

struct Fraction2 <: Fraction
    A::Vector{Float64}
    B::Vector{Float64}
    n::Int
end

struct Fraction3 <: Fraction
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    n::Int
end


function Fraction(endmembers::Int=2; n::Int=100)
    if endmembers==2
        Fraction(0, 1, n = n)
    elseif endmembers==3
        x = (0, 1/3)
        Fraction(x,x, n=n)
    else
        @error "Must provide either 2 or 3 endmembers."
    end
end

function Fraction(Amin::Number, Amax::Number; n::Int=100)
    @assert 0 ≤ Amin < Amax ≤ 1
    
    fA = collect(range(float(Amin), float(Amax), n))
    fB = 1.0 .- fA
    Fraction2(fA,fB,n)
end

function Fraction(A::Tuple{Number,Number}, B::Tuple{Number,Number}; n::Int=100)
        @assert 0 ≤ A[1] < A[2] ≤ 1
        @assert 0 ≤ B[1] < B[2] ≤ 1
        
        @assert A[2] + B[2] ≤ 1 


        fA = repeat(range(float(A[1]), float(A[2]), n), 1, n)
        fB = repeat(range(float(B[1]), float(B[2]), n), 1, n)'
        fC = @. 1 - fA + fB

        Fraction3(fA,fB,fC,n*n)
end

abstract type Model end

struct Model2 <: Model
    x::Vector{Float64}
    y::Vector{Float64}
end

struct Model3 <: Model
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

function Model(f::Fraction, s::T) where T<: System
    n = f.n
    c = countcomponents(s.A)

    if T <: System2 
        @assert c == countcomponents(s.B)
    elseif T<: System3
        @assert c == countcomponents(s.B) == countcomponents(s.C)
    end

    if c==4
        Model2(
            Vector{Float64}(undef,n),
            Vector{Float64}(undef,n),)
    elseif c==6
        Model3(
            Vector{Float64}(undef,n),
            Vector{Float64}(undef,n),
            Vector{Float64}(undef,n),)
    else
        @error "IsoMix does not support $(c/2) components"
    end
end


