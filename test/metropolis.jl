## initialguess, initialjump

testprior = Prior(Data(Norm(0,1),logNorm(2,1)), Data(Unf(2,6), Constant(8)), Data(Unconstrained(),Unconstrained()))

@test initialguess(testprior) isa System3{Component1}
@test initialjump(testprior) isa System3{Component1}

@test (initialguess(testprior).A.cx ≈ 7.38905609893065) & iszero(initialguess(testprior).A.x) & (initialguess(testprior).B ==  Component1(4.0, 8.0)) & (initialguess(testprior).C ==  Component1(0.5, 0.5))
@test (initialjump(testprior).A ==  Component1(1., 1.)) & (initialjump(testprior).B ==  Component1(1.,0.)) & (initialjump(testprior).C ==  Component1(0.1, 0.1))

## jump
testjump = jump(System(Component(1,2),Component(3,4)), System(Component(.1,.5),Component(.2,.4)), rng=StableRNG(1))
@test (testjump[1].A.x == 1.) & (testjump[1].A.cx == 2.) & (testjump[1].B.x != 3.)  & (testjump[1].B.cx == 4.)
@test testjump[2] == (:B, :x, 0.14332564801086908)


## update
@test update(System(Component(1,1),Component(1,1)), :B, :x, 12.) == System2{Component1}(Component1(1,1), Component1(12,1))
@test update(System(Component(1,1),Component(1,1), Component(1,1)), :B, :x, 12.) == System3{Component1}(Component1(1,1), Component1(12,1),Component1(1,1))

@test update(System(Component(1,1,1,1),Component(1,1,1,1)), :B, :x, 12.) == System2{Component2}(Component2(1,1,1,1), Component2(12,1,1,1))
@test update(System(Component(1,1,1,1),Component(1,1,1,1),Component(1,1,1,1)), :B, :x, 12.) == System3{Component2}(Component2(1,1,1,1), Component2(12,1,1,1),Component2(1,1,1,1))

@test update(System(Component(1,1,1,1,1,1),Component(1,1,1,1,1,1)), :B, :x, 12.) == System2{Component3}(Component3(1,1,1,1,1,1), Component3(12,1,1,1,1,1))
@test update(System(Component(1,1,1,1,1,1),Component(1,1,1,1,1,1),Component(1,1,1,1,1,1)), :B, :x, 12.) == System3{Component3}(Component3(1,1,1,1,1,1), Component3(12,1,1,1,1,1),Component3(1,1,1,1,1,1))