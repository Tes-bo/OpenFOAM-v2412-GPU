/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::incompressibleInterPhaseTransportModel

Description
    Transport model selection class for the incompressibleInterFoam family of
    solvers.

    By default the standard mixture transport modelling approach is used in
    which a single momentum stress model (laminar, non-Newtonian, LES or RAS) is
    constructed for the mixture.  However if the \c simulationType in
    constant/turbulenceProperties is set to \c twoPhaseTransport the alternative
    Euler-Euler two-phase transport modelling approach is used in which separate
    stress models (laminar, non-Newtonian, LES or RAS) are instantiated for each
    of the two phases allowing for different modeling for the phases.

SourceFiles
    incompressibleInterPhaseTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleInterPhaseTransportModel_H
#define incompressibleInterPhaseTransportModel_H


#include "turbulentTransportModel.H"
#include "VoFphaseTurbulentTransportModel.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class incompressibleInterPhaseTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class Mixture>
class incompressibleInterPhaseTransportModel
{
    // Private data

        //- Enum to select uniform or variable rho transport modelling
        enum rhoType
        {
            UNIFORM,
            VARIABLE
        };

        //- Rho type
        label rhoType_;

        //- Mixture volumetric flux
        const surfaceScalarField& phi_;

        //- Mixture mass flux
        const surfaceScalarField& rhoPhi_;

        //- Mixture incompressible turbulence with constant rho
        autoPtr<incompressible::turbulenceModel> incTurbulence_;

        //- Mixture incompressible turbulence with variable rho
        autoPtr<incompressible::phaseIncompressibleTurbulenceModel>
            rhoIncTurbulence_;


    // Private Member Functions

        //- No copy construct
        incompressibleInterPhaseTransportModel
        (
            const incompressibleInterPhaseTransportModel&
        ) = delete;

        //- No copy assignment
        void operator=(const incompressibleInterPhaseTransportModel&) = delete;


public:

     TypeName("incompressibleInterPhaseTransportModel");

    // Constructors

        //- Construct from components
        incompressibleInterPhaseTransportModel
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const surfaceScalarField& rhoPhi,
            const Mixture&
        );


    //- Destructor
    virtual ~incompressibleInterPhaseTransportModel() = default;


    // Member Functions

        //- Return the effective momentum stress divergence
        tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Correct the phase or mixture transport models
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "incompressibleInterPhaseTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
