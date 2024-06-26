# 2-endmember, 1 species
mtf, mts = Fraction(2,n=6), SysDraw(EmDraw(-3,2), EmDraw(3,4))
mtm = Model(mtf, mts)
mix!(mtm,mts,mtf)

@test mtm.x ≈ [3.0, 7/3, 1.5, 3/7, -1., -3.]
@test mtm.cx ≈ [4.0, 3.6, 3.2, 2.8, 2.4, 2.0]



# 2-endmember, 2 species
mtf, mts = Fraction(2,n=6), SysDraw(EmDraw(-3,2,-3,4), EmDraw(3,4,3,2))
mtm = mix(mts,mtf)

@test mtm.x ≈ [3.0, 7/3, 1.5, 3/7, -1., -3.]
@test mtm.y ≈ [3.0, 1.0, -3/7, -1.5, -7/3, -3.0]
@test mtm.cx ≈ [4.0, 3.6, 3.2, 2.8, 2.4, 2.0]
@test mtm.cy ≈ [2.0, 2.4000000000000004, 2.8, 3.2, 3.6, 4.0]


# 2-endmember, 3 species
mtf, mts = Fraction(2,n=6), SysDraw(EmDraw(-3,2,-3,4,-6,1), EmDraw(3,4,3,2,6,2))
mtm = mix(mts,mtf)

@test mtm.x ≈ [3.0, 7/3, 1.5, 3/7, -1., -3.]
@test mtm.y ≈ [3.0, 1.0, -3/7, -1.5, -7/3, -3.0]
@test mtm.z ≈ [6.0, 14/3, 3, 18/21, -2, -6.]
@test mtm.cx ≈ [4.0, 3.6, 3.2, 2.8, 2.4, 2.0]
@test mtm.cy ≈ [2.0, 2.4, 2.8, 3.2, 3.6, 4.0]
@test mtm.cz ≈ [2.0, 1.8, 1.6, 1.4, 1.2, 1.0]



# 3-endmember, 2 species
mtf, mts = Fraction(3,n=5), SysDraw(EmDraw(-3,2,-3,4), EmDraw(3,4,3,2), EmDraw(6,1,-6,1))
mtm = mix(mts,mtf)

@test mtm.x ≈ [6.0, 4.285714285714286, 3.6, 3.230769230769231, 3.0, 2.4, 2.25, 2.1818181818181817, 2.142857142857143, 0.0, 0.6666666666666666, 1.0, -1.7142857142857142, -0.6, -3.0]
@test mtm.y ≈ [-6.0, -2.4, 0.0, 1.7142857142857142, 3.0, -4.285714285714286, -2.25, -0.6666666666666666, 0.6, -3.6, -2.1818181818181817, -1.0, -3.230769230769231, -2.142857142857143, -3.0]
@test mtm.cx ≈ [1.0, 1.75, 2.5, 3.25, 4.0, 1.25, 2.0, 2.75, 3.5, 1.5, 2.25, 3.0, 1.75, 2.5, 2.0]
@test mtm.cy ≈ [1.0, 1.25, 1.5, 1.75, 2.0, 1.75, 2.0, 2.25, 2.5, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0]

# 3-endmember, 3 species
mtf, mts = Fraction(3,n=5), SysDraw(EmDraw(-3,2,-3,4,-6,1), EmDraw(3,4,3,2,6,2), EmDraw(6,1,-6,1, -2, 4))
mtm = mix(mts,mtf)

@test mtm.x ≈ [6.0, 4.285714285714286, 3.6, 3.230769230769231, 3.0, 2.4, 2.25, 2.1818181818181817, 2.142857142857143, 0.0, 0.6666666666666666, 1.0, -1.7142857142857142, -0.6, -3.0]
@test mtm.y ≈ [-6.0, -2.4, 0.0, 1.7142857142857142, 3.0, -4.285714285714286, -2.25, -0.6666666666666666, 0.6, -3.6, -2.1818181818181817, -1.0, -3.230769230769231, -2.142857142857143, -3.0]
@test mtm.z ≈ [-2.0, -0.8571428571428571, 0.6666666666666666, 2.8, 6.0, -2.3076923076923075, -0.9090909090909091, 1.1111111111111112, 4.285714285714286, -2.8, -1.0, 2.0, -3.7142857142857144, -1.2, -6.0]
@test mtm.cx ≈ [1.0, 1.75, 2.5, 3.25, 4.0, 1.25, 2.0, 2.75, 3.5, 1.5, 2.25, 3.0, 1.75, 2.5, 2.0]
@test mtm.cy ≈ [1.0, 1.25, 1.5, 1.75, 2.0, 1.75, 2.0, 2.25, 2.5, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0]
@test mtm.cz ≈ [4.0, 3.5, 3.0, 2.5, 2.0, 3.25, 2.75, 2.25, 1.75, 2.5, 2.0, 1.5, 1.75, 1.25, 1.0]


## Fraction calculators 

@test [fractions(SysDraw(EmDraw(-3,2,), EmDraw(3,4)), 7/3)...] ≈ [0.2,0.8]

@test isapprox([fractions(SysDraw(EmDraw(-3,2,-3,4), EmDraw(3,4,3,2), EmDraw(6,1,-6,1)), 3.6, 0.0)...], [0,.5,.5], atol = 1e-10)
