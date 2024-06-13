"""
    mix(s<:System, f<:Fraction)
"""
function mix(s::System, f::Fraction)
    m = Model(f,s)
    mix!(m,s,f)
    m
end


"""

    mix!(m <: Model, s <: System, f <: Fraction)

In-place version of [`mix`](@ref) that overwrites fields in `m`.

"""
function mix!(m::Model1, s::System2, f::Fraction2)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB = f.A[i], f.B[i]
        @fastmath Amass, Bmass = (fA * s.A.Cx, fB * s.B.Cx)
        @fastmath mixturemass =  Amass + Bmass
        @fastmath m.x[i] =   (Amass * s.A.x + Bmass * s.B.x ) / mixturemass
        m.Cx[i] = mixturemass
    end
    m
end

function mix!( m::Model2, s::System2, f::Fraction2)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB = f.A[i], f.B[i]
        @fastmath m.x[i] =   (fA * s.A.Cx * s.A.x + fB * s.B.Cx * s.B.x ) / (fA * s.A.Cx + fB * s.B.Cx)
        @fastmath m.y[i] = (fA * s.A.Cy * s.A.y + fB * s.B.Cy * s.B.y) / (fA * s.A.Cy + fB * s.B.Cy)
    end
    m
end

function mix!(m::Model3, s::System2, f::Fraction2)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB = f.A[i], f.B[i]
        @fastmath m.x[i] =   (fA * s.A.Cx * s.A.x + fB * s.B.Cx * s.B.x ) / (fA * s.A.Cx + fB * s.B.Cx)
        @fastmath m.y[i] = (fA * s.A.Cy * s.A.y + fB * s.B.Cy * s.B.y) / (fA * s.A.Cy + fB * s.B.Cy)
        @fastmath m.z[i] = (fA * s.A.Cz * s.A.z + fB * s.B.Cz * s.B.z) / (fA * s.A.Cz + fB * s.B.Cz)
    end
    m
end

function mix!(m::Model2, s::System3, f::Fraction3)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB, fC = f.A[i], f.B[i], f.C[i]
        @fastmath m.x[i] = (fA * s.A.Cx * s.A.x + fB * s.B.Cx * s.B.x + fC * s.C.Cx * s.C.x) / (fA * s.A.Cx + fB * s.B.Cx + fC * s.C.Cx)
        @fastmath m.y[i] = (fA * s.A.Cy * s.A.y + fB * s.B.Cy * s.B.y + fC * s.C.Cy * s.C.y) / (fA * s.A.Cy + fB * s.B.Cy + fC * s.C.Cy)
    end
    m
end

function mix!(m::Model3, s::System3, f::Fraction3)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB, fC = f.A[i], f.B[i], f.C[i]
        @fastmath m.x[i] = (fA * s.A.Cx * s.A.x + fB * s.B.Cx * s.B.x + fC * s.C.Cx * s.C.x) / (fA * s.A.Cx + fB * s.B.Cx + fC * s.C.Cx)
        @fastmath m.y[i] = (fA * s.A.Cy * s.A.y + fB * s.B.Cy * s.B.y + fC * s.C.Cy * s.C.y) / (fA * s.A.Cy + fB * s.B.Cy + fC * s.C.Cy)
        @fastmath m.z[i] = (fA * s.A.Cz * s.A.z + fB * s.B.Cz * s.B.z + fC * s.C.Cz * s.C.z) / (fA * s.A.Cz + fB * s.B.Cz + fC * s.C.Cz)
    end
    m
end