/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::radiation::wideBandDiffusiveRadiationMixedFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a wide-band, diffusive radiation
    condition, where the patch temperature is specified.

Usage
    \table
        Property     | Description             | Required    | Default value
        T            | temperature field name  | no          | T
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            wideBandDiffusiveRadiation;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::mixedFvPatchScalarField
    Foam::radiationCoupledBase

SourceFiles
    wideBandDiffusiveRadiationMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wideBandDiffusiveRadiationMixedFvPatchScalarField_H
#define wideBandDiffusiveRadiationMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
/*---------------------------------------------------------------------------*\
      Class wideBandDiffusiveRadiationMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wideBandDiffusiveRadiationMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("wideBandDiffusiveRadiation");


    // Constructors

        //- Construct from patch and internal field
        wideBandDiffusiveRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wideBandDiffusiveRadiationMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given GreyDiffusiveRadiationMixedFvPatchField
        //  onto a new patch
        wideBandDiffusiveRadiationMixedFvPatchScalarField
        (
            const wideBandDiffusiveRadiationMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wideBandDiffusiveRadiationMixedFvPatchScalarField
        (
            const wideBandDiffusiveRadiationMixedFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        wideBandDiffusiveRadiationMixedFvPatchScalarField
        (
            const wideBandDiffusiveRadiationMixedFvPatchScalarField&,
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

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
