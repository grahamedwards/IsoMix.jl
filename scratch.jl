## To-Do

using IsoMix, IsoMix.Examples


prior, measurements = Phillips2002()
p = IsoMix.initialguess(prior)

fraction = Fraction(3,n=101)

m = mix(p,fraction)

x=mixtropolis(prior,measurements,burninsteps=10_000, chainsteps=10_000)

using GLMakie

plot(x.ll)

####

s, f = SysDraw(EmDraw(-3,15,4,2,2,8), EmDraw(3,8,-2,22,5,12), EmDraw(-8,15,-10,16,12,3)), Fraction(3,n=1001)
m = mix(s, f)
