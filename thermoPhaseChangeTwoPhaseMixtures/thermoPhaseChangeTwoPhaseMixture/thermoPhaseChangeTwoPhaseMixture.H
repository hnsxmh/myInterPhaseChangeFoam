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
    Foam::thermoPhaseChangeTwoPhaseMixture

Description

SourceFiles
    thermoPhaseChangeTwoPhaseMixture.C
    phaseChangeModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef thermoPhaseChangeTwoPhaseMixture_H
#define thermoPhaseChangeTwoPhaseMixture_H

#include "thermoImmiscibleIncompressibleTwoPhaseMixture.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class thermoPhaseChangeTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class thermoPhaseChangeTwoPhaseMixture
:
    public thermoImmiscibleIncompressibleTwoPhaseMixture
{
protected:

    // Protected data

        dictionary thermoPhaseChangeTwoPhaseMixtureCoeffs_;

        // //- Saturation vapour pressure
        // dimensionedScalar pSat_;
        //- Latent heat
        dimensionedScalar lHeat_;        


public:

    //- Runtime type information
    TypeName("thermoPhaseChangeTwoPhaseMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermoPhaseChangeTwoPhaseMixture,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Constructors

        //- Construct from components
        thermoPhaseChangeTwoPhaseMixture
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        thermoPhaseChangeTwoPhaseMixture(const thermoPhaseChangeTwoPhaseMixture&);


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<thermoPhaseChangeTwoPhaseMixture> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~thermoPhaseChangeTwoPhaseMixture()
    {}


    // Member Functions

        //- Return const-access to the saturation vapour pressure
        // const volScalarField& pSat() const
        // {
        //     const volScalarField& pSat = alpha1_.db().lookupObject<volScalarField>("pSat");
        //     return pSat;
        // }

        //- Return const-access to the saturation vapour pressure
        const volScalarField& pSat() const
        {
            const volScalarField& pSat = DT_.db().lookupObject<volScalarField>("pSat");
            return pSat;
        }

        const dimensionedScalar& lHeat() const
        {
            return lHeat_;
        }

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField>> mDotP() const = 0;

        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        Pair<tmp<volScalarField>> vDotAlphal() const;

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients to multiply (p - pSat)
        Pair<tmp<volScalarField>> vDotP() const;

        //- Correct the phaseChange model
        virtual void correct() = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const thermoPhaseChangeTwoPhaseMixture&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
