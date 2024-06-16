export Phillips2002



Phillips2002() = (
    p = Prior(
        Data(Norm(-19.3,0.2), Norm(55,5),Norm(15.5,0.2), Norm(12,1)),
        Data(Norm(-16.6,0.2), Norm(51.5,5), Norm(7.9,0.2), Norm(14,1)),
        Data(Norm(-23.3,0.2), Norm(50.8,5), Norm(3.2,0.2), Norm(7.64,1))
    ),
    d = DataSet(x=Measurement([20.1],[0.2]),y=Measurement([7.6],[0.2]))
)