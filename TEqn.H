Pair<tmp<volScalarField>> vDotP = mixture.vDotP();
const volScalarField& vDotcP = vDotP[0]();
const volScalarField& vDotvP = vDotP[1]();
TSource = mixture.lHeat()/mixture.Cp()()*(vDotvP - vDotcP)*(mixture.pSat() - p_rgh);
{
    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
        + fvm::div(phi, T)
        - fvm::laplacian(DT, T)
        //when p_rgh is smaller than saturated vapor pressure, the liquid will evaporate and take away heat
        + TSource
        
    );
    TEqn.solve();
}
mixture.correct();
