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


@test -265 < loglikelihood(
    mix(SysDraw(EmDraw(-3,2), EmDraw(3,4)), Fraction(2,n=5)),
    Measurements(x=[-1 0.2; 3 0.2 ; 2 0.2], cx=[2 0.1; 1 0.1; NaN NaN])) < -264

testmodel = mix(SysDraw(EmDraw(-3,2,-3,4), EmDraw(3,4,3,2), EmDraw(6,1,-6,1)),Fraction(3,n=5))
testmsmt = Measurements([-1 0.2; 3 0.2 ; 2 0.2], [2 0.1; 1 0.1; NaN NaN])

@test -245 < loglikelihood(testmodel,testmsmt) < -244
@test -239 < loglikelihood(testmodel, Measurements([-1 0.2; 3 0.2 ; 2 0.2], [2 0.1; 1 0.1; 7 2])) < -237

@test -265 < loglikelihood(
    mix(SysDraw(EmDraw(-3,2), EmDraw(3,4)), Fraction(2,n=5)),
    Measurements(x=[-1 0.2; 3 0.2 ; 2 0.2], cx=[2 0.1; 1 0.1; NaN NaN])
    ) < -264

@test -253< loglikelihood(
    mix(SysDraw(EmDraw(-3,2,-3,4,2,1), EmDraw(3,4,3,2,-1,2), EmDraw(6,1,-6,1,4,6)),Fraction(3,n=5)),
    Measurements([-1 0.2; 3 0.2 ; 2 0.2], [2 0.1; 1 0.1; NaN NaN], [0.3 .2; NaN NaN; -0.4 .3])
    ) < -252

## nanll
@test IsoMix.nanll(.5, Norm(0,1)) == IsoMix.loglikelihood(.5, Norm(0,1))
@test iszero( IsoMix.nanll(NaN, Norm(0,1)) )

## top3

@test IsoMix.top3(-0.2, -Inf, -Inf, -Inf) == (-0.2, -Inf, -Inf)
@test IsoMix.top3(-0.2, -0.1, -Inf, -Inf) == (-0.1, -0.2, -Inf)
@test IsoMix.top3(-0.01, -0.1, -0.2, -Inf) == (-0.01, -0.1, -0.2)

@test IsoMix.top3(-8, -3, -4, -5) == (-3, -4, -5)
@test IsoMix.top3(-3.5, -3, -4, -3) == (-3, -3.5, -4,)
@test IsoMix.top3(-4.5, -3, -4, -5) == (-3, -4, -4.5)
@test IsoMix.top3(-2, -3, -4, -5) == (-2, -3, -4)

@test IsoMix.top3(-0.2, -0.1, -0.2, -Inf) == (-0.1, -0.2, -0.2)
@test IsoMix.top3(-0.1, -0.1, -0.2, -Inf) == (-0.1, -0.1, -0.2)
@test IsoMix.top3(0., -Inf, -Inf, -Inf) == (-Inf, -Inf, -Inf)