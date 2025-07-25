/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::temperaturePhaseChangeTwoPhaseMixture

Description

SourceFiles
    temperaturePhaseChangeTwoPhaseMixture.C
    newtemperaturePhaseChangeTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef temperaturePhaseChangeTwoPhaseMixture_H
#define temperaturePhaseChangeTwoPhaseMixture_H

#include "thermoIncompressibleTwoPhaseMixture.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class temperaturePhaseChangeTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class temperaturePhaseChangeTwoPhaseMixture
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the thermoIncompressibleTwoPhaseMixture
        const thermoIncompressibleTwoPhaseMixture& mixture_;

        //- Reference to fvMesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- No copy construct
        temperaturePhaseChangeTwoPhaseMixture
        (
            const temperaturePhaseChangeTwoPhaseMixture&
        ) = delete;

        //- No copy assignment
        void operator=(const temperaturePhaseChangeTwoPhaseMixture&) = delete;


public:

    //- Runtime type information
    TypeName("temperaturePhaseChangeTwoPhaseMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            temperaturePhaseChangeTwoPhaseMixture,
            components,
            (
                const thermoIncompressibleTwoPhaseMixture& mixture,
                const fvMesh& mesh
            ),
            (mixture, mesh)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<temperaturePhaseChangeTwoPhaseMixture> New
        (
            const thermoIncompressibleTwoPhaseMixture& mixture,
            const fvMesh& mesh
        );


    // Constructors

        //- Construct from components
        temperaturePhaseChangeTwoPhaseMixture
        (
            const thermoIncompressibleTwoPhaseMixture& mixture,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~temperaturePhaseChangeTwoPhaseMixture() = default;


    // Member Functions


        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        virtual Pair<tmp<volScalarField>> mDot() const = 0;

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (Tsat - T) for the condensation rate
        //  and a coefficient to multiply  (T - Tsat) for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotDeltaT() const = 0;

        //- Source for T equarion
        virtual tmp<fvScalarMatrix> TSource() const = 0;

        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> vDotAlphal() const;

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients
        virtual Pair<tmp<volScalarField>> vDot() const;

        //- Correct the phaseChange model
        virtual void correct() = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
