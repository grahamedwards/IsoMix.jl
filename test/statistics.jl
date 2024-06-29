## loglikelihood

@test loglikelihood(-1,Norm(0,1)) ≈ -0.5
@test loglikelihood(1,logNorm(2,1)) ≈ -2.
@test loglikelihood(2,Unf(0,1)) == -Inf
@test iszero(loglikelihood(0.5,Unf(0,1)))
@test iszero(loglikelihood(1,Constant(1)))
@test loglikelihood(1,Constant(1.2)) == -Inf
@test iszero(loglikelihood(rand(),Unconstrained()))

@test loglikelihood(EmDraw(1,1), Endmember(Norm(0,1),Norm(0,1))) ≈ -1.

@test -1.5 ≈ 
    loglikelihood(
        SysDraw(EmDraw(1,1), EmDraw(3,3)),
        Prior(Endmember(Norm(0,1), Norm(0,1)), Endmember(Norm(2,1), Unconstrained()))
)



testmodel = mix(SysDraw(EmDraw(-3,2,-3,4), EmDraw(3,4,3,2), EmDraw(6,1,-6,1)),Fraction(3,n=5))

testmsmt = Measurements([-1 0.2; 3 0.2 ; 2 0.2], [2 0.1; 1 0.1; NaN NaN])

@test -239 < loglikelihood(testmodel,testmsmt) < -238
@test -245 < loglikelihood(testmodel, Measurements([-1 0.2; 3 0.2 ; 2 0.2], [2 0.1; 1 0.1; 7 2])) < -244