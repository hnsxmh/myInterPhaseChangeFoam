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
Info << "T = [" << min(T).value() << ", " << max(T).value() << "]" << endl;
volScalarField CpConst
(
    IOobject
    (
        "TSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("", dimEnergy/dimMass/dimTemperature, 1.0)
);

dimensionedScalar kappaConst("kappa", dimensionSet(1,1,-3,-1,0,0,0), 1.0);

{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
        +   fvm::div(rhoPhi, T)
        // ==
        // 1.0 / CpConst * fvm::laplacian(kappaConst, T)
        // 1.0/mixture.Cp()()*fvm::laplacian(mixture.kappa()(), T)
        //when p_rgh is smaller than saturated vapor pressure, the liquid will evaporate and take away heat
        //- rho*TSource
    );
    TEqn.solve();
}
mixture.correct();
