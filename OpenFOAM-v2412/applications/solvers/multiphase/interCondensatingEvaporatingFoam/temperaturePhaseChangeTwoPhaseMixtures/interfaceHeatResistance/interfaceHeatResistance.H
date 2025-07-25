/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
    Copyright (C) 2020 Henning Scheufler
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
    Foam::phaseChangeTwoPhaseMixtures::interfaceHeatResistance

Description
    Interface Heat Resistance type of condensation/saturation model using
    spread source distribution following:

    References:
    \verbatim
        Hardt, S., Wondra, F. (2008).
        Evaporation model for interfacial flows based on a continuum-
        field representation of the source term
        Journal of Computational Physics 227 (2008), 5871-5895
    \endverbatim

SourceFiles
    interfaceHeatResistance.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceHeatResistance_H
#define interfaceHeatResistance_H

#include "temperaturePhaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace temperaturePhaseChangeTwoPhaseMixtures
{

/*---------------------------------------------------------------------------*\
                         Class interfaceHeatResistance
\*---------------------------------------------------------------------------*/

class interfaceHeatResistance
:
    public temperaturePhaseChangeTwoPhaseMixture
{
    // Private Data

        //- Heat transfer coefficient [1/s/K]
        dimensionedScalar R_;

        //- Interface area
        volScalarField interfaceArea_;

        //- Mass condensation source
        volScalarField mDotc_;

        //- Mass evaporation source
        volScalarField mDote_;

        //- Spread condensation mass source
        volScalarField mDotcSpread_;

        //- Spread evaporation mass source
        volScalarField mDoteSpread_;

        //- Spread for mass source
        scalar spread_;


    // Private Member Functions

        //- Update interface area
        void updateInterface();


public:

    //- Runtime type information
    TypeName("interfaceHeatResistance");


    // Constructors

        //- Construct from components
        interfaceHeatResistance
        (
            const thermoIncompressibleTwoPhaseMixture& mixture,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~interfaceHeatResistance() = default;


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        virtual Pair<tmp<volScalarField>> mDot() const;

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (Tsat - T) for the condensation rate
        //  and a coefficient to multiply  (T - Tsat) for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotDeltaT() const;

        //- Source for T equarion
        virtual tmp<fvScalarMatrix> TSource() const;

        //- Volumetric source for alpha (used by alphaEq)
        virtual Pair<tmp<volScalarField>> vDotAlphal() const;

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients (used by p_rghEq)
        virtual Pair<tmp<volScalarField>> vDot() const;

        //- Correct the interfaceHeatResistance phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace temperaturePhaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
