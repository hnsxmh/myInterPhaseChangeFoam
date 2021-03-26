Pair<tmp<volScalarField>> vDotP = mixture.vDotP();
const volScalarField& vDotcP = vDotP[0]();
const volScalarField& vDotvP = vDotP[1]();
massChange = (vDotvP - vDotcP)*(mixture.pSat() - p);
TSource = mixture.lHeat()/mixture.Cp()()*(vDotvP - vDotcP)*(mixture.pSat() - p);
Info << "Cp = [" << min(mixture.Cp()()).value() << ", " << max(mixture.Cp()()).value() << "]" << endl;
Info << "termal diffusivity = [" << min(mixture.DT()()).value() << ", " << max(mixture.DT()()).value() << "]" << endl;
Info << "mass change = [" << min(massChange).value() << ", " << max(massChange).value() << "]" << endl;
Info << "p_rgh = [" << min(p_rgh).value() << ", " << max(p_rgh).value() << "]" << endl;
Info << "p = [" << min(p).value() << ", " << max(p).value() << "]" << endl;

{
    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
        + fvm::div(phi, T)
        - fvm::laplacian(DT, T)
        //when p_rgh is smaller than saturated vapor pressure, the liquid will evaporate and take away heat
        - TSource
        
    );
    TEqn.solve();
}
mixture.correct();
