"""
    mix(s<:System, f<:Fraction)
"""
function mix(s::System, f::Fraction)
    m = Model(f,s)
    mix!(s,f,m)
    m
end


"""

    mix!(s <: System, f <: Fraction, m <: Model)

In-place version of [`mix`](@ref) that overwrites fields in `m`.

"""
function mix!(s::System2, f::Fraction2, m::Model1)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB = f.A[i], f.B[i]
        Amass, Bmass = (fA * s.A.Cx, fB * s.B.Cx)
        mixturemass =  Amass + Bmass
        m.x[i] =   (Amass * s.A.x + Bmass * s.B.x ) / mixturemass
        m.Cx[i] = mixturemass
    end
    m
end

function mix!(s::System2, f::Fraction2, m::Model2)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB = f.A[i], f.B[i]
        m.x[i] =   (fA * s.A.Cx * s.A.x + fB * s.B.Cx * s.B.x ) / (fA * s.A.Cx + fB * s.B.Cx)
        m.y[i] = (fA * s.A.Cy * s.A.y + fB * s.B.Cy * s.B.y) / (fA * s.A.Cy + fB * s.B.Cy)
    end
    m
end

function mix!(s::System2, f::Fraction2, m::Model3)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB = f.A[i], f.B[i]
        m.x[i] =   (fA * s.A.Cx * s.A.x + fB * s.B.Cx * s.B.x ) / (fA * s.A.Cx + fB * s.B.Cx)
        m.y[i] = (fA * s.A.Cy * s.A.y + fB * s.B.Cy * s.B.y) / (fA * s.A.Cy + fB * s.B.Cy)
        m.z[i] = (fA * s.A.Cz * s.A.z + fB * s.B.Cz * s.B.z) / (fA * s.A.Cz + fB * s.B.Cz)
    end
    m
end

function mix!(s::System3, f::Fraction3, m::Model2)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB, fC = f.A[i], f.B[i], f.C[i]
        m.x[i] = (fA * s.A.Cx * s.A.x + fB * s.B.Cx * s.B.x + fC * s.C.Cx * s.C.x) / (fA * s.A.Cx + fB * s.B.Cx + fC * s.C.Cx)
        m.y[i] = (fA * s.A.Cy * s.A.y + fB * s.B.Cy * s.B.y + fC * s.C.Cy * s.C.y) / (fA * s.A.Cy + fB * s.B.Cy + fC * s.C.Cy)
    end
    m
end

function mix!(s::System3, f::Fraction3, m::Model3)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB, fC = f.A[i], f.B[i], f.C[i]
        m.x[i] = (fA * s.A.Cx * s.A.x + fB * s.B.Cx * s.B.x + fC * s.C.Cx * s.C.x) / (fA * s.A.Cx + fB * s.B.Cx + fC * s.C.Cx)
        m.y[i] = (fA * s.A.Cy * s.A.y + fB * s.B.Cy * s.B.y + fC * s.C.Cy * s.C.y) / (fA * s.A.Cy + fB * s.B.Cy + fC * s.C.Cy)
        m.z[i] = (fA * s.A.Cz * s.A.z + fB * s.B.Cz * s.B.z + fC * s.C.Cz * s.C.z) / (fA * s.A.Cz + fB * s.B.Cz + fC * s.C.Cz)
    end
    m
end