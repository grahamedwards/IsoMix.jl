"""
    mix(s<:System, f<:Fraction)

Calculate mixing between components of system `s` given fractional mixtures `f`, returned in a `Model` instance.

see also: [`mix`](@ref), [`System`](@ref), [`Fraction`](@ref), [`Model`](@ref)

---

### Input configurations:

    mix(s::System2{Component1}, f::Fraction2) -> Model1
Two-endmember, 1 species.

    mix(s::System2{Component2}, f::Fraction2) -> Model2
Two-endmember, 2 species.

    mix(s::System2{Component3}, f::Fraction2) -> Model3
Two-endmember, 3 species.

    mix(s::System3{Component2}, f::Fraction3) -> Model2
Three-endmember, 2 species.

    mix(s::System3{Component3}, f::Fraction3) -> Model3
Three-endmember, 3 species.


"""
function mix(s::System, f::Fraction)
    m = Model(f,s)
    mix!(m,s,f)
    m
end



"""

    mix!(m <: Model, s <: System, f <: Fraction)

In-place version of [`mix`](@ref) that overwrites fields in `m`.

see also: [`System`](@ref), [`Fraction`](@ref), [`Model`](@ref)

"""
function mix!(m::Model1, s::System2, f::Fraction2)
    @inbounds @simd ivdep for i = 1:f.n
        fxA, fxB = f.A[i] * s.A.cx, f.B[i] * s.B.cx

        xmix = fxA+fxB

        m.x[i] = (fxA * s.A.x + fxB * s.B.x) / xmix
        m.cx[i] = xmix
    end
    m
end

function mix!( m::Model2, s::System2, f::Fraction2)
    @inbounds @simd ivdep for i = 1:f.n
        fxA, fyA = f.A[i] .* (s.A.cx, s.A.cy)
        fxB, fyB = f.B[i] .* (s.B.cx, s.B.cy)

        xmix, ymix = fxA+fxB, fyA+fyB

        m.x[i] = (fxA * s.A.x + fxB * s.B.x) / xmix
        m.y[i] = (fyA * s.A.y + fyB * s.B.y) / ymix

        m.cx[i], m.cy[i] = xmix, ymix
    end
    m
end

function mix!(m::Model3, s::System2, f::Fraction2)
    @inbounds @simd ivdep for i = 1:f.n

        fxA, fyA, fzA = f.A[i] .* (s.A.cx, s.A.cy, s.A.cz)
        fxB, fyB, fzB = f.B[i] .* (s.B.cx, s.B.cy, s.B.cz)

        xmix, ymix, zmix = fxA+fxB, fyA+fyB, fzA+fzB

        m.x[i] = (fxA * s.A.x + fxB * s.B.x) / xmix
        m.y[i] = (fyA * s.A.y + fyB * s.B.y) / ymix
        m.z[i] = (fzA * s.A.z + fzB * s.B.z) / zmix

        m.cx[i], m.cy[i], m.cz[i] = xmix, ymix, zmix
    end
    m
end

function mix!(m::Model2, s::System3, f::Fraction3)
    @inbounds @simd ivdep for i = 1:f.n

        fxA, fyA = f.A[i] .* (s.A.cx, s.A.cy)
        fxB, fyB = f.B[i] .* (s.B.cx, s.B.cy)
        fxC, fyC = f.C[i] .* (s.C.cx, s.C.cy)

        xmix, ymix = fxA+fxB+fxC, fyA+fyB+fyC

        m.x[i] = (fxA * s.A.x + fxB * s.B.x + fxC * s.C.x) / xmix
        m.y[i] = (fyA * s.A.y + fyB * s.B.y + fyC * s.C.y) / ymix

        m.cx[i], m.cy[i] = xmix, ymix
    end
    m
end

function mix!(m::Model3, s::System3, f::Fraction3)
    @inbounds @simd ivdep for i = 1:f.n

        fxA, fyA, fzA = f.A[i] .* (s.A.cx, s.A.cy, s.A.cz)
        fxB, fyB, fzB = f.B[i] .* (s.B.cx, s.B.cy, s.B.cz)
        fxC, fyC, fzC= f.C[i] .* (s.C.cx, s.C.cy, s.C.cz)

        xmix, ymix, zmix = fxA+fxB+fxC, fyA+fyB+fyC, fzA+fzB+fzC

        m.x[i] = (fxA * s.A.x + fxB * s.B.x + fxC * s.C.x) / xmix
        m.y[i] = (fyA * s.A.y + fyB * s.B.y + fyC * s.C.y) / ymix
        m.z[i] = (fzA * s.A.z + fzB * s.B.z + fzC * s.C.z) / zmix
        
        m.cx[i], m.cy[i], m.cz[i] = xmix, ymix, zmix 
    end
    m
end



"""

    fractions(s<:System2{Component1}, x)

Calculate the fractional contributions of each the two `Component`s in `s`, given a composition of `x`.

---

    fractions(s<:System3{Component2}, x, y)

Calculate the fractional contributions of each the three `Component`s in `s`, given compositions of `x` and `y`.

"""
function fractions(s::System2{Component1}, x::Float64) 
    out = StaticArrays.SA[(s.A.x-x)*s.A.cx (s.B.x-x)*s.B.cx; 1 1] \ StaticArrays.SA[0,1]
    (out.x,out.y)
end

function fractions(s::System3{Component2}, x::T, y::T) where T<: Float64
    out = StaticArrays.SA[(s.A.x-x)*s.A.cx (s.B.x-x)*s.B.cx (s.C.x-x)*s.C.cx ; (s.A.y-y)*s.A.cx (s.B.y-y)*s.B.cy (s.C.y-y)*s.C.cy ; 1 1 1] \ StaticArrays.SA[0 , 0 , 1]
    (out.x,out.y,out.z)
end