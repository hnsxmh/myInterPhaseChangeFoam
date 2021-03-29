{
const dimensionedScalar t30_11("30.11", dimensionSet(0,0,0,1,0,0,0), 30.11);
const dimensionedScalar t273_15("273.15", dimensionSet(0,0,0,1,0,0,0), 273.15);
const dimensionedScalar t1("1", dimensionSet(0,0,0,1,0,0,0), 1);
const dimensionedScalar p610_94("610.94", dimensionSet(1,-1,-2,0,0,0,0), 610.94);
// dimensionSet( [kg], [m], [s], [K], [kg*mol], [A], [cd]), [kg/(m*S^2)]=[Pa]
// August-Roche-Magnus formula
// pSat = p610_94 * exp( 17.625*(T-t273_15) / max(t1, T-t30_11) );
pSat = dimensionedScalar("pSatconst", dimPressure, 2300);
Info << "pSat = [" << min(pSat).value() << ", " << max(pSat).value() << "]" << endl;
//max(1,...) is included to avoid problems with devision by 0
}