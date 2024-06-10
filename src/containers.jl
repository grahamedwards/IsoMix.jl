

struct Component
    x::Float64
    y::Float64
    z::Float64
end

Component(x::Number,y::Number) = Component(x,y,NaN)

function countcomponents(c::Component)
    f = fieldnames(Component)
    n = 0

    @inbounds for i in f
        n += ifelse(isnan(getfield(c,i)),0,1)
    end
    n
end


abstract type System end

struct System2 <: System
    A::Component
    cA::Component
    B::Component
    cB::Component
end

struct System3 <: System
    A::Component
    cA::Component
    B::Component
    cB::Component
    C::Component
    cC::Component
end

System(A::T, B::T; cA::T=T(1,1),cB::T=T(1,1)) where T <: Component= System2(A,cA,B,cB)
System(A::T, B::T, C::T; cA::T=T(1,1),cB::T=T(1,1), cC::T = T(1,1)) where T <: Component = System3(A,cA,B,cB,C,cC)


abstract type Prior end 

struct Norm <: Prior
    mu::Float64
    sig::Float64
end

struct logNorm <: Prior
    mu::Float64
    sig::Float64
end

struct Unf <: Prior
    a::Float64
    b::Float64
end

struct Constant <: Prior
    x::Float64
end

abstract type Data end

struct Data2 <: Data
    x::Prior
    cx::Prior
    y::Prior
    cy::Prior
end

struct Data3 <: Data
    x::Prior
    cx::Prior
    y::Prior
    cy::Prior
    z::Prior
    cz::Prior
end


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

        Fraction3(fA,fB,fC,n)
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

    @assert c == countcomponents(s.cA) ==
    countcomponents(s.B) ==  countcomponents(s.cB)

    T <: System3 && @assert c == countcomponents(s.C) ==  countcomponents(s.cC)

    if c==2
        Model2(
            Vector{Float64}(undef,n),
            Vector{Float64}(undef,n),)
    elseif c==3
        Model3(
            Vector{Float64}(undef,n),
            Vector{Float64}(undef,n),
            Vector{Float64}(undef,n),)
    else
        @error "IsoMix does not support $c components"
    end
end


