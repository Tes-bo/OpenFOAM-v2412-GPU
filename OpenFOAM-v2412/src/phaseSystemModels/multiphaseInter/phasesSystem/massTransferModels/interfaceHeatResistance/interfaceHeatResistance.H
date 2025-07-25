/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Henning Scheufler
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::meltingEvaporationModels::interfaceHeatResistance

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


Usage

    Example usage:
    \verbatim
    massTransferModel
    (
        (liquid to gas)
        {
            type                interfaceHeatResistance;
            R                   2e6;
            Tactivate           373;
        }
    );
    \endverbatim

    where:
    \table
        Property    | Description             | Required    | Default value
        R           | Heat transfer coefficient | yes
        includeVolChange    | Volumen change  | no          | yes
        isoAlpha    | iso-alpha for interface | no          | 0.5
        Tactivate   | Saturation temperature  | yes
        species     | Specie name on the other phase | no   | none
        spread      | Cells to spread the source for pEq | no   | 3
    \endtable

SourceFiles
    interfaceHeatResistance.C

\*---------------------------------------------------------------------------*/

#ifndef meltingEvaporationModels_interfaceHeatResistance_H
#define meltingEvaporationModels_interfaceHeatResistance_H

#include "InterfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{

// Forward Declarations
class phasePair;

namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                   Class interfaceHeatResistance Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class interfaceHeatResistance
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private Data

        //- Heat transfer coefficient [1/s/K]
        dimensionedScalar R_;

        //- Activation temperature
        const dimensionedScalar Tactivate_;

        //- Interface area
        volScalarField interfaceArea_;

        //- Mass source
        volScalarField mDotc_;

        //- Spread mass source
        volScalarField mDotcSpread_;

        //- Heat transfer coefficient
        volScalarField htc_;

        //- Interface Iso-value
        scalar isoAlpha_;

        //- Spread for mass source
        scalar spread_;


    // Private Member Functions

        //- Update interface
        void updateInterface(const volScalarField& T);

public:

    //- Runtime type information
    TypeName("interfaceHeatResistance");


    // Constructors

        //- Construct from components
        interfaceHeatResistance
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interfaceHeatResistance() = default;


    // Member Functions

        //- Explicit total mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            const volScalarField& field
        );

        //- Implicit mass transfer coefficient
        virtual tmp<volScalarField> KSp
        (
            label modelVariable,
            const volScalarField& field
        );

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> KSu
        (
            label modelVariable,
            const volScalarField& field
        );

        //- Return Tactivate
        virtual const dimensionedScalar& Tactivate() const noexcept
        {
            return Tactivate_;
        }

        //- Add/subtract alpha*div(U) as a source term
        //- for alpha, substituting div(U) = mDot(1/rho1 - 1/rho2)
        virtual bool includeDivU() const noexcept
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interfaceHeatResistance.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
