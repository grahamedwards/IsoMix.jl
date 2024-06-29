import IsoMix
export Phillips2002



Phillips2002() = (
    prior = IsoMix.Prior(
        IsoMix.Data(IsoMix.Norm(-19.3,0.2), IsoMix.Norm(55,5),IsoMix.Norm(15.5,0.2), IsoMix.Norm(12,1)),
        IsoMix.Data(IsoMix.Norm(-16.6,0.2), IsoMix.Norm(51.5,5), IsoMix.Norm(7.9,0.2), IsoMix.Norm(14,1)),
        IsoMix.Data(IsoMix.Norm(-23.3,0.2), IsoMix.Norm(50.8,5), IsoMix.Norm(3.2,0.2), IsoMix.Norm(7.64,1))
    ),
    measurements = IsoMix.Measurements([20.1 0.2], [7.6 0.2])
)



