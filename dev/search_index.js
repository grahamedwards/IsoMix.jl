var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = IsoMix","category":"page"},{"location":"#IsoMix","page":"Home","title":"IsoMix","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for IsoMix.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [IsoMix]","category":"page"},{"location":"#IsoMix.Calculations","page":"Home","title":"IsoMix.Calculations","text":"Calculations <: IsoMixType\n\nAbstract supertype for types used specifically in chemical and isotopic model simulations. \n\nsee also: Fraction, Model\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Constant","page":"Home","title":"IsoMix.Constant","text":"Constant(x) <: Datum\n\nA Datum instance with a discrete value x.\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Data","page":"Home","title":"IsoMix.Data","text":"Data <: IsoMixType\n\nAbstract supertype for types that contain measured data, represented by Datum instances. Prior instances contain 2-3 Endmember instances, and both Endmember orMeasurementsinstances containDatum` instances.\n\nsee also: Datum, Endmember, Prior, Measurements\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Datum","page":"Home","title":"IsoMix.Datum","text":"Datum <: Data\n\nAbstract supertype for Datum instances.\n\nsee also: Constant, Unf, Norm, logNorm, Unconstrained\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.DistributionDraws","page":"Home","title":"IsoMix.DistributionDraws","text":"DistributionDraws <: IsoMixType\n\nAbstract supertype for types that contain values drawn from Data distributions of Endmember and Prior instances. SysDraw instances contain 2-3 EmDraw instances, which contain discrete draws from a corresponding Endmember instance.\n\nsee also: EmDraw, SysDraw\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.EmDraw","page":"Home","title":"IsoMix.EmDraw","text":"EmDraw <: DistributionDraws\n\nShort for EndmemberDraw. Abstract supertype for EmDraw_ instances, where the _ indicates dimensionality. Represents the composition of an endmember within a natural system. \n\nsee also: SysDraw, EmDraw1, EmDraw2, EmDraw3\n\n\n\nEmDraw(x, cx) -> EmDraw1\n\nEmDraw(x, cx, y, cy) -> EmDraw2\n\nEmDraw(x, cx, y, cy, z, cz) -> EmDraw3\n\nThe constructor function accepts values as a list of arguments or as keyword assignments, e.g.: EmDraw(x=1, cx=2)\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.EmDraw1","page":"Home","title":"IsoMix.EmDraw1","text":"EmDraw1 <: EmDraw\n\n1-dimensional EmDraw instance.\n\nFields \nx Isotopic composition of x\ncx Concentration of x\n\nsee also: EmDraw\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.EmDraw2","page":"Home","title":"IsoMix.EmDraw2","text":"EmDraw2 <: EmDraw\n\n2-dimensional EmDraw instance.\n\nFields \nx Isotopic composition of x\ncx Concentration of x\ny Isotopic composition of y\ncy Concentration of y\n\nsee also: EmDraw\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.EmDraw3","page":"Home","title":"IsoMix.EmDraw3","text":"EmDraw3 <: EmDraw\n\n3-dimensional EmDraw instance.\n\nFields \nx Isotopic composition of x\ncx Concentration of x\ny Isotopic composition of y\ncy Concentration of y\nz Isotopic composition of z\ncz Concentration of z\n\nsee also: EmDraw\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Endmember","page":"Home","title":"IsoMix.Endmember","text":"Endmember <: Data\n\nAbstract supertype for Endmember_ instances, where the _ indicates dimensionality.\n\nsee also: Endmember1, Endmember2, Endmember3\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Endmember1","page":"Home","title":"IsoMix.Endmember1","text":"Endmember1 <: Endmember\n\n1-dimensional Endmember instance.\n\nFields \nx Isotopic composition of x\ncx Concentration of x\n\nsee also: Endmember, Datum\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Endmember2","page":"Home","title":"IsoMix.Endmember2","text":"Endmember2 <: Endmember\n\n2-dimensional Endmember instance.\n\nFields \nx Isotopic composition of x\ncx Concentration of x\ny Isotopic composition of y\ncy Concentration of y\n\nsee also: Endmember, Datum\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Endmember3","page":"Home","title":"IsoMix.Endmember3","text":"Endmember3 <: Endmember\n\n3-dimensional Endmember instance.\n\nFields \nx Isotopic composition of x\ncx Concentration of x\ny Isotopic composition of y\ncy Concentration of y\nz Isotopic composition of z\ncz Concentration of z\n\nsee also: Endmember, Datum\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Fraction","page":"Home","title":"IsoMix.Fraction","text":"Fraction <: Calculations\n\nAbstract supertype for Fraction_ instances that contain vectors of fractional combinations of _ components/endmembers.\n\nsee also: Fraction2, Fraction3\n\n\n\nFraction(Amin, Amax; n=101)\n\nReturns a Fraction2 with n linearly spaced fractional combinations of components A and B that each sum to 1, given minimum and maximum A fractions of Amin and Amax, respectively. Sincen gives the linear fractions of each component, this constructor returns n-length vectors of each component. Both Amin and Amax must fall within [0,1].\n\nExample\n\njulia> Fraction(0,1,n=3)\nFraction2([0.0, 0.5, 1.0], [1.0, 0.5, 0.0], 3)\n\n\n\nFraction(Amin, Amax, Bmin, Bmax; n=101)\n\nReturns a Fraction3 with linearly spaced fractional combinations of components A, B, and C that each sum to 1, given respective minima and maxima of A (Amin and Amax) and B (Bmin and Bmax). Since n gives the linear fractions of each component, this constructor returns a n(n÷2+1)-length vector for each component. Amin,Amax, Bmin, and Bmax must fall within [0,1].\n\nExample\n\njulia> f=Fraction(0,1,0,1,n=3)\nFraction3([0.0, 0.0, 0.0, 0.5, 0.5, 1.0], [0.0, 0.5, 1.0, 0.0, 0.5, 0.0], [1.0, 0.5, 0.0, 0.5, 0.0, 0.0], 6)\n\n\n\nFraction(dim, n=101)\n\nGiven a provided number of components/dimensions dim (2 or 3), returns a corresponding Fraction instance with component concentrations of 0 to 1. \n\nExample\n\njulia> Fraction(2,n=3)\nFraction2([0.0, 0.5, 1.0], [1.0, 0.5, 0.0], 3)\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Fraction2","page":"Home","title":"IsoMix.Fraction2","text":"Fraction2 <: Fraction\n\n2-dimensional Fraction instance.\n\nFields \nA Endmember/component A\nB Endmember/component B\n\nsee also: Fraction\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Fraction3","page":"Home","title":"IsoMix.Fraction3","text":"Fraction3 <: Fraction\n\n3-dimensional Fraction instance.\n\nFields \nA Endmember/component A\nB Endmember/component B\nC Endmember/component C\n\nsee also: Fraction\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.IsoMixType","page":"Home","title":"IsoMix.IsoMixType","text":"IsoMixType\n\nSupertype containing all custom types and structs in the IsoMix package:\n\nDirect subtypes: EmDraw, Endmember, Measurements, Datum, Fraction, Measurements, Model,Prior, SysDraw \n\n\n\nIsoMixType\n├─ Calculations\n│  ├─ Fraction\n│  │  ├─ Fraction2\n│  │  └─ Fraction3\n│  └─ Model\n│     ├─ Model1\n│     ├─ Model2\n│     └─ Model3\n├─ Data\n│  ├─ Datum\n│  │  ├─ Constant\n│  │  ├─ Norm\n│  │  ├─ Unconstrained\n│  │  ├─ Unf\n│  │  └─ logNorm\n│  ├─ Endmember\n│  │  ├─ Endmember1\n│  │  ├─ Endmember2\n│  │  └─ Endmember3\n│  ├─ Measurements\n│  │  ├─ Measurements1\n│  │  ├─ Measurements2\n│  │  └─ Measurements3\n│  └─ Prior\n│     ├─ Prior2\n│     └─ Prior3\n└─ DistributionDraws\n   ├─ EmDraw\n   │  ├─ EmDraw1\n   │  ├─ EmDraw2\n   │  └─ EmDraw3\n   └─ SysDraw\n      ├─ SysDraw2\n      └─ SysDraw3\n\n\n\nGenerated with \n\nusing AbstractTrees,IsoMix; AbstractTrees.children(d) = subtypes(d); print_tree(IsoMixType)\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Measurements","page":"Home","title":"IsoMix.Measurements","text":"Measurements <: Data\n\nAbstract supertype for Measurements_ instances, where _ indicates dimensionality. Represents a dataset of measured (normally distributed) data as vectors of Norm instances. Absent Measurement fields are represented as vectors of Unconstrained and do not affect loglikelihood calculations. \n\nsee also: Measurements, Measurements1, Measurements2, Measurements3\n\nConstruction\n\nThe constructor function Measurements(...) takes nx2 matrices, where n represents the number of measured samples in the dataset, and n must be the same for all provided matrices. Each matrix represents a series of measurements (e.g. isotopic composition or abundance measurement) for a given species. The table below defines all possible fields of a Measurements and which Measurements subtypes contain them (denoted with a ✓):\n\nFields  Measurements1 Measurements2 Measurements3\nx Isotopic composition of x ✓ ✓ ✓\ncx Concentration of x ✓ ✓ ✓\n    \ny Isotopic composition of y  ✓ ✓\ncy Concentration of y  ✓ ✓\n    \nz Isotopic composition of z   ✓\ncz Concentration of z   ✓\n\nNote that each row among the Measurement fields corresponds to a specific sample, so provided matrices must consistently reflect the same sample in each row. Measurements with missing values must be replaced with NaN. \n\nInputing data into Measurements()\n\nThe constructor function accepts a list of Measurements instances:\n\nWithout concentration measurements:\n\nMeasurements(x, y) -> Measurements2 \nMeasurements(x, y, z) -> Measurements3\n\nWith concentration measurements:\n\nMeasurements(x=..., cx=...) -> Measurements1\nMeasurements(x, cx, y, cy) -> Measurements2\nMeasurements(x, cx, y, cy, z, cz) -> Measurements3\n\nFor all other datasets configurations, declare measurements with keywords:\n\nMeasurements(; x, cx, y, cy, z, cz)\n\n...which determines the appropriate subtype and all undeclared fields return vectors of Unconstrained instances.\n\nExamples\n\n\njulia> Measurements(x = [1 2; 3 4], cx = [5 6; 7 8])\nMeasurements1(Datum[Norm(1.0, 2.0), Norm(3.0, 4.0)], Datum[Norm(5.0, 6.0), Norm(7.0, 8.0)])\n\njulia> m2d = Measurements([1 2; 3 4], [5 6; 7 NaN]); \n\njulia> m2d.x\n2-element Vector{Datum}:\n Norm(1.0, 2.0)\n Norm(3.0, 4.0)\n\njulia> m2d.cx\n2-element Vector{Datum}:\n Unconstrained()\n Unconstrained()\n\njulia> m2d.y\n2-element Vector{Datum}:\n Norm(5.0, 6.0)\n Norm(7.0, NaN)\n\njulia> m2d = Measurements([1 2; 3 4], [5 6; 7 8; 9 10])\nERROR: AssertionError: x & y must have the same number of entries (rows). Use NaN for missing values.\n\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Measurements1","page":"Home","title":"IsoMix.Measurements1","text":"Measurements1 <: Measurements\n\nMeasurements instance reflecting isotope and abundance measurements of 1 species: x.\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Measurements2","page":"Home","title":"IsoMix.Measurements2","text":"Measurements2 <: Measurements\n\nMeasurements instance reflecting isotope and abundance measurements of 2 species: x, y.\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Measurements3","page":"Home","title":"IsoMix.Measurements3","text":"Measurements3 <: Measurements\n\nMeasurements instance reflecting isotope and abundance measurements of 3 species: x, y, z.\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Model","page":"Home","title":"IsoMix.Model","text":"Model <: Calculations\n\nAbstract supertype for Model_ instances, where the _ indicates dimensionality.\n\nsee also: Model1, Model2, Model3\n\n\n\nModel(f::Fraction, s::SysDraw)\n\nGenerate a Model instance \n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Model1","page":"Home","title":"IsoMix.Model1","text":"Model1 <: Model\n\n1-dimensional Model instance: a single species tracking isotopic composition and concentration.\n\nFields \nx Isotopic composition of x\ncx Concentration of x\n\nsee also: Model\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Model2","page":"Home","title":"IsoMix.Model2","text":"Model2 <: Model\n\n2-dimensional Model instance: two species tracking isotopic composition and concentration.\n\nFields \nx Isotopic composition of x\ncx Concentration of x\ny Isotopic composition of y\ncy Concentration of y\n\nsee also: Model\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Model3","page":"Home","title":"IsoMix.Model3","text":"Model3 <: Model\n\n3-dimensional Model instance: three species tracking isotopic composition and concentration.\n\nFields \nx Isotopic composition of x\ncx Concentration of x\ny Isotopic composition of y\ncy Concentration of y\nz Isotopic composition of z\ncz Concentration of z\n\nsee also: Model\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Norm","page":"Home","title":"IsoMix.Norm","text":"Norm(m, s) <: Datum\n\nNormally distributed Datum with mean m and standard deviation s.\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Prior","page":"Home","title":"IsoMix.Prior","text":"Prior <: Data\n\nAbstract supertype for Prior_ instances, where the _ indicates dimensionality. Represents a suite of [Endmember] data.\n\nsee also: Endmember, Prior2, Prior3\n\n\n\nPrior(A::Endmember, B::Endmember)\n\nReturns a Prior2\n\nPrior(A::Endmember, B::Endmember, C::Endmember)\n\nReturns a Prior3\n\nNote: each Endmember must be of the same subtype (i.e. consistent number of components within each endmember of a system).\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Prior2","page":"Home","title":"IsoMix.Prior2","text":"Prior2 <: Prior\n\n2-dimensional Prior instance.\n\nFields \nA Endmember/component A\nB Endmember/component B\n\nsee also: Prior\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Prior3","page":"Home","title":"IsoMix.Prior3","text":"Prior3 <: Prior\n\n3-dimensional Prior instance.\n\nFields \nA Endmember/component A\nB Endmember/component B\nC Endmember/component C\n\nsee also: Prior\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.SysDraw","page":"Home","title":"IsoMix.SysDraw","text":"SysDraw <: Data\n\nShort for SystemDraw. Abstract supertype for SysDraw_ instances, where the _ indicates dimensionality. Represents a model system composition charaterized by its EmDraw fields.\n\nsee also: EmDraw SysDraw2, SysDraw3\n\n\n\nThe constructor function accepts values as a list of arguments:\n\nSysDraw(A::EmDraw, B::EmDraw) -> SysDraw2\n\nSysDraw(A::EmDraw, B::EmDraw, C::EmDraw) -> SysDraw3\n\n... or as keyword assignments, e.g.:      SysDraw(A=..., B=...)\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.SysDraw2","page":"Home","title":"IsoMix.SysDraw2","text":"SysDraw2 <: SysDraw\n\n2-dimensional SysDraw instance.\n\nFields \nA Endmember/component A\nB Endmember/component B\n\nsee also: SysDraw\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.SysDraw3","page":"Home","title":"IsoMix.SysDraw3","text":"SysDraw3 <: SysDraw\n\n3-dimensional SysDraw instance.\n\nFields \nA Endmember/component A\nB Endmember/component B\nC Endmember/component C\n\nsee also: SysDraw\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Unconstrained","page":"Home","title":"IsoMix.Unconstrained","text":"Unconstrained <: Datum\n\nA Datum instance for an unconstrained variable, functionally similar to Unf(-∞,∞) or an unmeasured Measurements datum. Typically, using this value is inadvisable.\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.Unf","page":"Home","title":"IsoMix.Unf","text":"Unf(a, b) <: Datum\n\nUniformly distributed Datum with minimum a and maximum b, inclusive. I.e., for any x in Unf(a b), x in ab.\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.logNorm","page":"Home","title":"IsoMix.logNorm","text":"logNorm(lm, ls) <: Datum\n\nLog-normally distributed Datum with log-mean lm and log-space standard deviation ls.\n\n\n\n\n\n","category":"type"},{"location":"#IsoMix.countcomponents-Tuple{C} where C<:EmDraw","page":"Home","title":"IsoMix.countcomponents","text":"IsoMix.countcomponents(c::EmDraw)\n\nCount the number of elemental dimensions (1-3) of c.\n\nExample\n\njulia> IsoMix.countcomponents(EmDraw(1,2,3,4,5,6))\n3\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.extractcomponents-Tuple{EmDraw1}","page":"Home","title":"IsoMix.extractcomponents","text":"IsoMix.extractcomponents(c)\n\nReturn all values from EmDraw instance c.\n\nsee also: IsoMix.extractsystem\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.extractfields-Union{Tuple{S}, Tuple{C}} where {C<:EmDraw, S<:SysDraw{C}}","page":"Home","title":"IsoMix.extractfields","text":"IsoMix.extractfields(s)\n\nReturn all SysDraw field names (e.g. :A, :B) and all EmDraw subfield names (e.g. :x, :cx) as single Symbols defining each, e.g. :Ax, :Bcz.\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.extractsystem-Tuple{SysDraw3}","page":"Home","title":"IsoMix.extractsystem","text":"IsoMix.extractsystem(s)\n\nReturns all EmDraw field values within the SysDraw fields of s.\n\nsee also: IsoMix.extractcomponents\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.fractions-Tuple{SysDraw2{EmDraw1}, Float64}","page":"Home","title":"IsoMix.fractions","text":"fractions(s<:SysDraw2{EmDraw1}, x)\n\nCalculate the fractional contributions of each the two EmDraws in s, given a composition of x.\n\n\n\nfractions(s<:SysDraw3{EmDraw2}, x, y)\n\nCalculate the fractional contributions of each the three EmDraws in s, given compositions of x and y.\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.initialguess-Union{Tuple{P}, Tuple{D}} where {D<:Endmember, P<:Prior{D}}","page":"Home","title":"IsoMix.initialguess","text":"IsoMix.initialguess(p<:Prior)\n\nReturns a SysDraw instance corresponding to provided Prior instance p, assigning initial guesses of each component composition using the central tendancy of each Datum: mean of Norm, log-mean of logNorm, midpoint of Unf, value of Constant, and an arbitrary value of 0.5 for Unconstrained.\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.initialjump-Union{Tuple{P}, Tuple{D}} where {D<:Endmember, P<:Prior{D}}","page":"Home","title":"IsoMix.initialjump","text":"IsoMix.initialjump(p<:Prior)\n\nReturns a SysDraw instance corresponding to an initial jumping distribution σ for each component's compositions, given provided Prior instance p. Jump σ are assigned for each Datum as follows, the standar deviation (σ) of Norm, the log-σ of logNorm, ¼ the range of Unf, a value of 0 for Constant, and an arbitrary value of 0.1 for Unconstrained.\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.jump-Union{Tuple{S}, Tuple{C}, Tuple{S, S}} where {C<:EmDraw, S<:SysDraw{C}}","page":"Home","title":"IsoMix.jump","text":"IsoMix.jump(s<:SysDraw, j<:SysDraw; rng)\n\nGiven a guess of SysDraw values s and corresponding jumping distribution scales in j, randomly perturbs one component of s and returns the new guess as well as a tuple containing the jumped EmDraw/endmember (:A, :B, :C), the jumped component composition (e.g. :x, :y, :cx), and the jump value.\n\nExample\n\njulia> jumpedsystem, t = jump(SysDraw(EmDraw(1,2),EmDraw(3,4)), SysDraw(EmDraw(.1,.5),EmDraw(.2,.4)), rng=IsoMix.Random.Xoshiro(1)); \n\njulia> jumpedsystem\nSysDraw2{EmDraw1}(EmDraw1(1.0, 2.349413341845734), EmDraw1(3.0, 4.0))\n\njulia> t\n(:A, :cx, 0.34941334184573425)\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.loglikelihood-Tuple{Number, Norm}","page":"Home","title":"IsoMix.loglikelihood","text":"loglikelihood(x::Number, D <: Datum)\n\nCalculate the relative¹ log-likelihood that x is drawn from a distribution D, which may be normal (Norm), lognormal (logNorm), or uniform (Unf). Note that log-space constants are dropped from the calculations for efficiency (metropolis only compares log-ratios)\n\n\n\nloglikelihood(c<:EmDraw, d<:Endmember)\n\nCalculate the loglikelihood that the species composition/concentrations in c were drawn from the corresponding distributions in d.\n\n\n\nloglikelihood(s<:SysDraw , p<:Prior)\n\nCalculate the loglikelihood that the components in s were drawn from the corresponding Endmember in p.\n\n\n\nloglikelihood(model<:Model , measurements<:Measurements)\n\nCalculate the likelihood (ℓ) that the (highest likelihood) model compositions were drawn from the distributions in measurements. Uses Polyester.@batch-based multithreading for faster runtimes. \n\nNote that rather than calculating the ℓ of every model composition relative to the corresponding distribution in 'measurements, the mean of the three highest ℓ among all instances in amodelare summed among measurements. This preventsmixtropolisfrom converging on small parameter spaces that only just fit around themeasurements`, while also limiting aliasing problems.\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.mix!-Tuple{Model1, SysDraw2, Fraction2}","page":"Home","title":"IsoMix.mix!","text":"mix!(m <: Model, s <: SysDraw, f <: Fraction)\n\nIn-place version of mix that overwrites fields in m.\n\nsee also: SysDraw, Fraction, Model\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.mix-Tuple{SysDraw, Fraction}","page":"Home","title":"IsoMix.mix","text":"mix(s<:SysDraw, f<:Fraction)\n\nCalculate mixing between components of system s given fractional mixtures f, returned in a Model instance.\n\nsee also: mix, SysDraw, Fraction, Model\n\n\n\nInput configurations:\n\nmix(s::SysDraw2{EmDraw1}, f::Fraction2) -> Model1\n\nTwo-endmember, 1 species.\n\nmix(s::SysDraw2{EmDraw2}, f::Fraction2) -> Model2\n\nTwo-endmember, 2 species.\n\nmix(s::SysDraw2{EmDraw3}, f::Fraction2) -> Model3\n\nTwo-endmember, 3 species.\n\nmix(s::SysDraw3{EmDraw2}, f::Fraction3) -> Model2\n\nThree-endmember, 2 species.\n\nmix(s::SysDraw3{EmDraw3}, f::Fraction3) -> Model3\n\nThree-endmember, 3 species.\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.nanll-Tuple{Float64, Norm}","page":"Home","title":"IsoMix.nanll","text":"IsoMix.nanll(x,D<:Union{Norm,Unconstrained})\n\nSame as loglikelihood, but returns 0.0 instead of NaN. Returns 1 for Unconstrained priors.\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.reassigncomponents-Union{Tuple{T}, Tuple{Symbol, T, T, T}} where T<:EmDraw","page":"Home","title":"IsoMix.reassigncomponents","text":"IsoMix.reassigncomponents(si::Symbol, X, A, B)\nIsoMix.reassigncomponents(si::Symbol, X, A, B,C)\n\nReturns a SysDraw containing EmDraws A and B (and C), replacing the field denoted by si with X\n\nExamples\n\njulia> IsoMix.reassigncomponents(:A, EmDraw(2,2), EmDraw(1,1), EmDraw(1,1))\n\nSysDraw2{EmDraw1}(EmDraw1(2.0, 2.0), EmDraw1(1.0, 1.0))\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.stopwatch-Tuple{Integer, Integer, Number}","page":"Home","title":"IsoMix.stopwatch","text":"IsoMix.stopwatch(i, n, t)\n\nConvenience function for [porewatermetropolis] that returns a String reporting the progress at step i for total steps n with start time t (in s since the epoch).\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.top3-NTuple{4, Number}","page":"Home","title":"IsoMix.top3","text":"IsoMix.top3(y, x1, x2, x3)\n\nReturns the three highest inputs in decreasing order, as long as x1, x2, and x3 are provided in decreasing order. This is optimized for speed and a very specific use case, so entering x values out of order will result in inaccurate ranking.\n\n\n\n\n\n","category":"method"},{"location":"#IsoMix.update","page":"Home","title":"IsoMix.update","text":"update(s<:SysDraw, si::Symbol, ci::Symbol, v)\n\nUpdate the EmDraw field ci of SysDraw field si of s with the value v.\n\nExample\n\njulia> update(SysDraw(EmDraw(1,2),EmDraw(3,4)), :B, :x, 12.)\nSysDraw2{EmDraw1}(EmDraw1(1.0, 2.0), EmDraw1(12.0, 4.0))\n\n\n\n\n\n","category":"function"}]
}
