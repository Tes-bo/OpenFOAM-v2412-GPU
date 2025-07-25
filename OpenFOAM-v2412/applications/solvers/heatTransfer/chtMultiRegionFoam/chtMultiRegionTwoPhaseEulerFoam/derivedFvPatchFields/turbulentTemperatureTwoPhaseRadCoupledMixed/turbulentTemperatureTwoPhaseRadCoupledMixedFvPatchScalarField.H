/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd
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
    Foam::compressible::
        turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature and radiation heat transfer
    to be used for in multiregion cases with two phase Euler system

Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field    | no | T
        qrNbr      | name of the radiative flux in the nbr region | no | none
        qr         | name of the radiative flux in this region | no | none
        region     | region to which this BC belongs | yes
        otherPhase | name of the vapour phase in the fluid region | yes
        kappaMethod  | inherited from temperatureCoupledBase | inherited |
        kappa        | inherited from temperatureCoupledBase | inherited |
    \endtable

    Example of the boundary condition specification on the fluid region:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureTwoPhaseRadCoupledMixed;
        Tnbr            T;
        qrNbr           none;
        qr              none;
        kappaMethod     phaseSystem;
        region          fluid;
        otherPhase      gas;
        value           uniform 300;
    }
    \endverbatim

    Example of the boundary condition specification on the solid region:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureTwoPhaseRadCoupledMixed;
        Tnbr            T.liquid;
        qrNbr           none;
        qr              none;
        kappaMethod     solidThermo;
        region          solid;
        otherPhase      gas;
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.


SourceFiles
    turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField_H
#define turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:

     // Public enumerations

        //- Type of supplied Kappa
        enum KMethodType
        {
            mtSolidThermo,
            mtLookup,
            mtPhaseSystem
        };

     // Data types

        //- Enumeration listing the region
        enum regionType
        {
            solid,
            fluid
        };

private:

    // Private data

        //- Heat source type names
        static const Enum<regionType> regionTypeNames_;

        //- Kappa method types
        static const Enum<KMethodType> KMethodTypeNames_;

        //- Heat source type
        regionType regionType_;

        //- How to get K
        const KMethodType method_;

        //- Name of thermal conductivity field (if looked up from database)
        const word kappaName_;

        //- name of the other phase (vapor/liquid phase)
        word otherPhaseName_;

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Name of the radiative heat flux in the neighbour region
        const word qrNbrName_;

        //- Name of the radiative heat flux in local region
        const word qrName_;


    // Private members

        //- Given patch temperature calculate corresponding K field
        tmp<scalarField> kappa(const scalarField& Tp) const;


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureTwoPhaseRadCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
        (
            const
            turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
        (
            const turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
