## loglikelihood

@test loglikelihood(-1,Norm(0,1)) ≈ -0.5
@test loglikelihood(1,logNorm(2,1)) ≈ -2.
@test loglikelihood(2,Unf(0,1)) == -Inf
@test iszero(loglikelihood(0.5,Unf(0,1)))
@test iszero(loglikelihood(1,Constant(1)))
@test loglikelihood(1,Constant(1.2)) == -Inf
@test iszero(loglikelihood(rand(),Unconstrained()))

@test loglikelihood(Component(1,1), Data(Norm(0,1),Norm(0,1))) ≈ -1.

testsys = 
testprr = 

@test -1.5 ≈ 
    loglikelihood(
        System(Component(1,1), Component(3,3)),
        Prior(Data(Norm(0,1), Norm(0,1)), Data(Norm(2,1), Unconstrained()))
)

@test iszero(loglikelihood([1.,2.,3.,4.], Measurements()))

@test loglikelihood([1.,2.,3.,4.], Measurements([1,4],[.5,2])) ≈ -29.75

@test loglikelihood([1.,2.,3.,4.], Measurements([1,NaN],[.5,2])) ≈ loglikelihood([1.,2.,3.,4.], Measurements([1],[0.5]))


testmodel = mix(System(Component(-3,2,-3,4), Component(3,4,3,2), Component(6,1,-6,1)),Fraction(3,n=5))

testdataset = DataSet(
    Measurements([-1, 3,2],[0.2,0.2,0.2]),
    Measurements([2,1,2.5],[0.1,0.1,0.1]),
    Measurements([-2,-3,0],[0.3,0.3,0.3]),
    Measurements()
)

@test loglikelihood(testmodel,testdataset) ≈ -9256.200833693123