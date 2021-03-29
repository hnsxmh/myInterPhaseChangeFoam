/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::thermoIncompressibleTwoPhaseMixture

Description
    A two-phase incompressible transportModel

SourceFiles
    thermoIncompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef thermoIncompressibleTwoPhaseMixture_H
#define thermoIncompressibleTwoPhaseMixture_H

#include "transportModel.H"
#include "viscosityModel.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class thermoIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class thermoIncompressibleTwoPhaseMixture
:
    public IOdictionary,
    public transportModel,
    public twoPhaseMixture
{
protected:

    // Protected data

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;
        // ///////////////////////
        dimensionedScalar Cp1_;
        dimensionedScalar Cp2_;

        dimensionedScalar kappa1_;
        dimensionedScalar kappa2_;

        ////////////////////////
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        volScalarField nu_;

        //- Laminar thermal diffusivity [m2/s]/////[kg/m/s]
        volScalarField DT_;

    // Protected Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();
        void calcDT();


public:

    TypeName("thermoIncompressibleTwoPhaseMixture");


    // Constructors

        //- Construct from components
        thermoIncompressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~thermoIncompressibleTwoPhaseMixture()
    {}


    // Member Functions

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        /////////////////////////////////////////////////////////////

        const dimensionedScalar& Cp1() const
        {
            return Cp1_;
        };

        const dimensionedScalar& Cp2() const
        {
            return Cp2_;
        };

        const dimensionedScalar& kappa1() const
        {
            return kappa1_;
        };

        const dimensionedScalar& kappa2() const
        {
            return kappa2_;
        };

        // //- Return the thermal diffusivity
        const volScalarField& DT() const
        {
            return DT_;
        }

        ///////////////////////////////////////////////////////////////        

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- Return the Cp
        virtual tmp<volScalarField> Cp() const;

        //- Return the kappa
        virtual tmp<volScalarField> kappa() const;

        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcNu();
            calcDT();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
