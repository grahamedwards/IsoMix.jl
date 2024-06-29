## IsoMix.initialguess, IsoMix.initialjump

testprior = Prior(Endmember(Norm(0,1),logNorm(2,1)), Endmember(Unf(2,6), Constant(8)), Endmember(Unconstrained(),Unconstrained()))

@test IsoMix.initialguess(testprior) isa SysDraw3{EmDraw1}
@test IsoMix.initialjump(testprior) isa SysDraw3{EmDraw1}

@test (IsoMix.initialguess(testprior).A.cx ≈ 7.38905609893065) & iszero(IsoMix.initialguess(testprior).A.x) & (IsoMix.initialguess(testprior).B ==  EmDraw1(4.0, 8.0)) & (IsoMix.initialguess(testprior).C ==  EmDraw1(0.5, 0.5))
@test (IsoMix.initialjump(testprior).A ==  EmDraw1(1., 1.)) & (IsoMix.initialjump(testprior).B ==  EmDraw1(1.,0.)) & (IsoMix.initialjump(testprior).C ==  EmDraw1(0.1, 0.1))

## jump
testjump = IsoMix.jump(SysDraw(EmDraw(1,2),EmDraw(3,4)), SysDraw(EmDraw(.1,.5),EmDraw(.2,.4)), rng=StableRNG(1))
@test (testjump[1].A.x == 1.) & (testjump[1].A.cx == 2.) & (testjump[1].B.x != 3.)  & (testjump[1].B.cx == 4.)
@test testjump[2] == (:B, :x, 0.14332564801086908)


## update
@test update(SysDraw(EmDraw(1,1),EmDraw(1,1)), :B, :x, 12.) == SysDraw2{EmDraw1}(EmDraw1(1,1), EmDraw1(12,1))
@test update(SysDraw(EmDraw(1,1),EmDraw(1,1), EmDraw(1,1)), :B, :x, 12.) == SysDraw3{EmDraw1}(EmDraw1(1,1), EmDraw1(12,1),EmDraw1(1,1))

@test update(SysDraw(EmDraw(1,1,1,1),EmDraw(1,1,1,1)), :B, :x, 12.) == SysDraw2{EmDraw2}(EmDraw2(1,1,1,1), EmDraw2(12,1,1,1))
@test update(SysDraw(EmDraw(1,1,1,1),EmDraw(1,1,1,1),EmDraw(1,1,1,1)), :B, :x, 12.) == SysDraw3{EmDraw2}(EmDraw2(1,1,1,1), EmDraw2(12,1,1,1),EmDraw2(1,1,1,1))

@test update(SysDraw(EmDraw(1,1,1,1,1,1),EmDraw(1,1,1,1,1,1)), :B, :x, 12.) == SysDraw2{EmDraw3}(EmDraw3(1,1,1,1,1,1), EmDraw3(12,1,1,1,1,1))
@test update(SysDraw(EmDraw(1,1,1,1,1,1),EmDraw(1,1,1,1,1,1),EmDraw(1,1,1,1,1,1)), :B, :x, 12.) == SysDraw3{EmDraw3}(EmDraw3(1,1,1,1,1,1), EmDraw3(12,1,1,1,1,1),EmDraw3(1,1,1,1,1,1))


@test IsoMix.extractsystem(SysDraw3(EmDraw(1,2,3,4,5,6), EmDraw(1,2,3,4,5,6), EmDraw(1,2,3,4,5,6))) == (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

@test IsoMix.extractsystem(SysDraw2(EmDraw(1,2,3,4,5,6), EmDraw(1,2,3,4,5,6))) == (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
@test IsoMix.extractsystem(SysDraw2(EmDraw(1,2,3,4), EmDraw(1,2,3,4))) == (1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0)
@test IsoMix.extractsystem(SysDraw2(EmDraw(1,2,), EmDraw(1,2))) == (1.0, 2.0, 1.0, 2.0)