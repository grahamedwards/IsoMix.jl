## First test type structure

@test EmDraw <: IsoMixType
@test (EmDraw1 <: EmDraw) & (EmDraw2 <: EmDraw) & (EmDraw3 <: EmDraw)

@test Endmember <: IsoMixType
@test (Endmember1 <: Endmember) & (Endmember2 <: Endmember) & (Endmember3 <: Endmember)

@test Measurements <: IsoMixType
@test (Measurements1 <: Measurements) & (Measurements2 <: Measurements) & (Measurements3 <: Measurements)

@test Datum <: IsoMixType
@test prod((Constant, Norm, logNorm, Unf, Unconstrained) .<: Datum)

@test Fraction <: IsoMixType
@test (Fraction2 <: Fraction) & (Fraction3 <: Fraction)

@test Model <: IsoMixType
@test (Model1 <: Model) & (Model2 <: Model) & (Model3 <: Model)

@test Measurements <: IsoMixType

@test Prior <: IsoMixType
@test (Prior2 <: Prior) & (Prior3 <: Prior)

@test SysDraw <: IsoMixType
@test (SysDraw2 <: SysDraw) & (SysDraw3 <: SysDraw)

## Test constructors 

# EmDraw

@test EmDraw(1,1) === EmDraw1(1,1)
@test EmDraw(1,1,1,1) === EmDraw2(1,1,1,1)
@test EmDraw(1,1,1,1,1,1) === EmDraw3(1,1,1,1,1,1)

@test EmDraw(x=1,cx=1) === EmDraw1(1,1)
@test EmDraw(y=1,cy=1,x=1,cx=1) === EmDraw2(1,1,1,1)
@test EmDraw(x=1,y=1,z=1,cx=1,cy=1,cz=1) === EmDraw3(1,1,1,1,1,1)

@test (IsoMix.countcomponents(EmDraw(1,1)), IsoMix.countcomponents(EmDraw(1,1,1,1))) === (1,2)

# SysDraw

@test SysDraw(EmDraw(1,1),EmDraw(2,2)) == SysDraw2(EmDraw(1,1),EmDraw(2,2))

@test SysDraw(EmDraw(1,1),EmDraw(2,2),EmDraw(3,3)) == SysDraw3(EmDraw(1,1),EmDraw(2,2),EmDraw(3,3))

# Datum

x = Norm(1,2); @test (x.m,x.s) === (1.,2.)
x = logNorm(1,2); @test (x.lm,x.ls) === (1.,2.)
x = Unf(1,2); @test (x.a,x.b) === (1.,2.)
x = Constant(1); @test x.x === 1.

# Endmember

@test Endmember(Constant(1),Norm(1,1)) === Endmember1(Constant(1),Norm(1,1))
@test Endmember(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1)) === Endmember2(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1))
@test Endmember(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1),Unf(1,1),Unf(1,1)) === Endmember3(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1),Unf(1,1),Unf(1,1))

@test Endmember(x=Constant(1),cx=Norm(1,1)) === Endmember1(Constant(1),Norm(1,1))
@test Endmember(x=Constant(1),cx=Norm(1,1),y=logNorm(1,1),cy=Unf(1,1)) === Endmember2(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1))
@test Endmember(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1),Unf(1,1),Unf(1,1)) === Endmember3(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1),Unf(1,1),Unf(1,1))

# Prior

@test Prior(Endmember1(Constant(1),Norm(1,1)),Endmember1(Norm(1,1),Constant(1))) == Prior2(Endmember1(Constant(1),Norm(1,1)),Endmember1(Norm(1,1),Constant(1)))

@test Prior(Endmember1(Constant(1),Norm(1,1)),Endmember(Norm(1,1),Constant(1)),Endmember1(Norm(1,1),Unf(0,1))) == Prior3(Endmember1(Constant(1),Norm(1,1)),Endmember1(Norm(1,1),Constant(1)),Endmember1(Norm(1,1),Unf(0,1)))

@test try Prior(Endmember1(Constant(1),Norm(1,1)),Endmember2(Norm(1,1),Constant(1), Norm(1,1),Unf(0,1))) catch e true end

# Fraction 

testfraction = Fraction(2,n=3)
@test typeof(testfraction) <: Fraction2

testfracbool = 1
testfracbool *= isapprox(testfraction.A, [0.0, 0.5, 1.0])
testfracbool *= isapprox(testfraction.B, [1.0, 0.5, 0.0])
@test testfracbool * testfraction.n == 3
@test μ(testfraction.A .+ testfraction.B) ≈ 1.0

testfraction = Fraction(3, n=3)
@test typeof(testfraction) <: Fraction3

testfracbool *= isapprox(testfraction.A, [0.0, 0.0, 0.0, 0.5, 0.5, 1.0])
testfracbool *= isapprox(testfraction.B, [0.0, 0.5, 1.0, 0.0, 0.5, 0.0])
testfracbool *= isapprox(testfraction.C, [1.0, 0.5, 0.0, 0.5, 0.0, 0.0])
@test testfracbool * testfraction.n == 6
@test μ(@. testfraction.A + testfraction.B + testfraction.C) ≈ 1.0

@test isapprox(Fraction(.2,.6,n=3).B, [0.8, 0.6, 0.4]) 

@test isapprox(Fraction(.2,.6,.3,.7,n=3).B, [0.3, 0.5, 0.7, 0.3, 0.5, 0.3]) 

# Model

@test length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1),EmDraw(2,2,))).cx) == length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1),EmDraw(2,2))).x) ==3 

@test length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1,1,1),EmDraw(2,2,2,2))).y) == length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1,1,1),EmDraw(2,2,2,2))).x) == length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1,1,1),EmDraw(2,2,2,2))).cy) == length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1,1,1),EmDraw(2,2,2,2))).cx) ==3 

@test length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1,1,1,1,1),EmDraw(2,2,2,2,2,2))).z) == length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1,1,1,1,1),EmDraw(2,2,2,2,2,2))).y) == length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1,1,1,1,1),EmDraw(2,2,2,2,2,2))).x) == length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1,1,1,1,1),EmDraw(2,2,2,2,2,2))).cz) == length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1,1,1,1,1),EmDraw(2,2,2,2,2,2))).cy) == length(Model(Fraction(2,n=3), SysDraw(EmDraw(1,1,1,1,1,1),EmDraw(2,2,2,2,2,2))).cx) == 3 

# Measurements

@test Measurements(x = [1 2; 3 4], cx = [5 6; 7 8]) isa Measurements1
@test Measurements(x = [1 2; 3 4], cx = [5 6; 7 8]).x == Datum[Norm(1.0, 2.0), Norm(3.0, 4.0)]

@test Measurements([1 2; 3 4], [5 6; 7 8]) isa Measurements2
@test Measurements([1 2; 3 4], [5 6; 7 8]).cx == Datum[Unconstrained(), Unconstrained()]
@test Measurements([1 2; 3 4], [5 6; 7 8]).y == Datum[Norm(5.0, 6.0), Norm(7.0, 8.0)]

@test Measurements([1 2; 3 4], [5 6; 7 8], [9 10; 11 12]) isa Measurements3
@test Measurements([1 2; 3 4], [5 6; 7 8], [9 10; 11 12]).cz == Datum[Unconstrained(), Unconstrained()]
@test Measurements([1 2; 3 4], [5 6; 7 8], [9 10; 11 12]).z == Datum[Norm(9.0, 10.0), Norm(11.0, 12.0)]

@test Measurements([1 2; 3 4], [5 6; 7 8], [9 10; 11 12], [13 14; 15 16]) isa Measurements2
@test Measurements([1 2; 3 4], [5 6; 7 8], [9 10; 11 12], [13 14; 15 16], [17 18; 19 20], [21 22; 23 24]) isa Measurements3