


function mix(s::System; f::Fraction=Fraction())
    m = Model(f,s)
    mix!(x,f,m)
    m
end





function mix!(s::System2, f::Fraction2, m::Model2)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB = f.A[i], f.B[i]
        m.x[i] =   (fA * s.cA.x * s.A.x + fB * s.cB.x * s.B.x ) / (fA * s.cA.x + fB * s.cB.x)
        m.y[i] = (fA * s.cA.y * s.A.y + fB * s.cB.y * s.B.y) / (fA * s.cA.y + fB * s.cB.y)
    end
    m
end

function mix!(s::System2, f::Fraction2, m::Model3)
    @inbounds @simd ivdep for i = 1:f.n
        fA, fB = f.A[i], f.B[i]
        m.x[i] =   (fA * s.cA.x * s.A.x + fB * s.cB.x * s.B.x ) / (fA * s.cA.x + fB * s.cB.x)
        m.y[i] = (fA * s.cA.y * s.A.y + fB * s.cB.y * s.B.y) / (fA * s.cA.y + fB * s.cB.y)
        m.z[i] = (fA * s.cA.z * s.A.z + fB * s.cB.z * s.B.z) / (fA * s.cA.z + fB * s.cB.z)
    end
    m
end
