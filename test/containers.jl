## First test type structure

@test Component <: IsoMixType
@test (Component1 <: Component) & (Component2 <: Component) & (Component3 <: Component)

@test Data <: IsoMixType
@test (Data1 <: Data) & (Data2 <: Data) & (Data3 <: Data)

@test Datum <: IsoMixType
@test (Constant <: Datum) & (Norm <: Datum) & (logNorm <: Datum) & (Unf <: Datum)

@test Fraction <: IsoMixType
@test (Fraction2 <: Fraction) & (Fraction3 <: Fraction)

@test Model <: IsoMixType
@test (Model1 <: Model) & (Model2 <: Model) & (Model3 <: Model)

@test Prior <: IsoMixType
@test (Prior2 <: Prior) & (Prior3 <: Prior)

@test System <: IsoMixType
@test (System2 <: System) & (System3 <: System)

## Test constructors 

# Component

@test Component(1,1) === Component1(1,1)
@test Component(1,1,1,1) === Component2(1,1,1,1)
@test Component(1,1,1,1,1,1) === Component3(1,1,1,1,1,1)

@test Component(x=1,Cx=1) === Component1(1,1)
@test Component(y=1,Cy=1,x=1,Cx=1) === Component2(1,1,1,1)
@test Component(x=1,y=1,z=1,Cx=1,Cy=1,Cz=1) === Component3(1,1,1,1,1,1)

@test (IsoMix.countcomponents(Component(1,1)), IsoMix.countcomponents(Component(1,1,1,1))) === (1,2)

# System

@test System(Component(1,1),Component(2,2)) == System2(Component(1,1),Component(2,2))

@test System(Component(1,1),Component(2,2),Component(3,3)) == System3(Component(1,1),Component(2,2),Component(3,3))

# Datum

x = Norm(1,2); @test (x.m,x.s) === (1.,2.)
x = logNorm(1,2); @test (x.lm,x.ls) === (1.,2.)
x = Unf(1,2); @test (x.a,x.b) === (1.,2.)
x = Constant(1); @test x.x === 1.

# Data

@test Data(Constant(1),Norm(1,1)) === Data1(Constant(1),Norm(1,1))
@test Data(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1)) === Data2(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1))
@test Data(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1),Unf(1,1),Unf(1,1)) === Data3(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1),Unf(1,1),Unf(1,1))

@test Data(x=Constant(1),Cx=Norm(1,1)) === Data1(Constant(1),Norm(1,1))
@test Data(x=Constant(1),Cx=Norm(1,1),y=logNorm(1,1),Cy=Unf(1,1)) === Data2(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1))
@test Data(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1),Unf(1,1),Unf(1,1)) === Data3(Constant(1),Norm(1,1),logNorm(1,1),Unf(1,1),Unf(1,1),Unf(1,1))

# Prior

@test Prior(Data1(Constant(1),Norm(1,1)),Data1(Norm(1,1),Constant(1))) == Prior2(Data1(Constant(1),Norm(1,1)),Data1(Norm(1,1),Constant(1)))

@test Prior(Data1(Constant(1),Norm(1,1)),Data(Norm(1,1),Constant(1)),Data1(Norm(1,1),Unf(0,1))) == Prior3(Data1(Constant(1),Norm(1,1)),Data1(Norm(1,1),Constant(1)),Data1(Norm(1,1),Unf(0,1)))

@test try Prior(Data1(Constant(1),Norm(1,1)),Data2(Norm(1,1),Constant(1), Norm(1,1),Unf(0,1))) catch e true end

# Fraction 

testfraction = Fraction(2,n=3)
@test typeof(testfraction) <: Fraction2

testfracbool = 1
testfracbool *= prod(testfraction.A .≈ [0.0, 0.5, 1.0])
testfracbool *= prod(testfraction.B .≈ [1.0, 0.5, 0.0])
@test testfracbool * testfraction.n == 3
@test μ(testfraction.A .+ testfraction.B) ≈ 1.0

testfraction = Fraction(3, n=3)
@test typeof(testfraction) <: Fraction3

testfracbool *= prod(testfraction.A .≈ [0.0, 0.0, 0.0, 0.5, 0.5, 1.0])
testfracbool *= prod(testfraction.B .≈ [0.0, 0.5, 1.0, 0.0, 0.5, 0.0])
testfracbool *= prod(testfraction.C .≈ [1.0, 0.5, 0.0, 0.5, 0.0, 0.0])
@test testfracbool * testfraction.n == 6
@test μ(@. testfraction.A + testfraction.B + testfraction.C) ≈ 1.0

@test prod(Fraction(.2,.6,n=3).B .≈ (0.8, 0.6, 0.4)) 

@test prod(Fraction(.2,.6,.3,.7,n=3).B .≈ (0.3, 0.5, 0.7, 0.3, 0.5, 0.3)) 

# Model


@test length(Model(Fraction(2,n=3), System(Component(1,1),Component(2,2,))).Cx) == length(Model(Fraction(2,n=3), System(Component(1,1),Component(2,2))).x) ==3 

@test length(Model(Fraction(2,n=3), System(Component(1,1,1,1),Component(2,2,2,2))).y) == length(Model(Fraction(2,n=3), System(Component(1,1,1,1),Component(2,2,2,2))).x) ==3 

@test length(Model(Fraction(2,n=3), System(Component(1,1,1,1,1,1),Component(2,2,2,2,2,2))).z) == length(Model(Fraction(2,n=3), System(Component(1,1,1,1,1,1),Component(2,2,2,2,2,2))).y) == length(Model(Fraction(2,n=3), System(Component(1,1,1,1,1,1),Component(2,2,2,2,2,2))).x) == 3 

