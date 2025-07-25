/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::compressible::vibrationShellFvPatchScalarField

Group
    grpVibrationBoundaryConditions

Description
    This boundary condition provides a coupled acoustic pressure condition
    between a primary region (3D mesh) and a vibration shell model (2D mesh).

Usage
    Example of the boundary condition specification:
    \verbatim
    <masterPatchName>
    {
        // Mandatory entries
        type                vibrationShell;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                      | Type  | Reqd | Deflt
      type     | Type name: vibrationShell        | word  | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link mixedFvPatchField.H \endlink
      - \link vibrationShellModel.H \endlink

SourceFiles
    vibrationShellFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef vibrationShellFvPatchScalarField_H
#define vibrationShellFvPatchScalarField_H

#include "vibrationShellModel.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class vibrationShellFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class vibrationShellFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
    // Typedefs

        //- The finite-area region model
        typedef regionModels::vibrationShellModel baffleType;


    // Private Data

        //- The vibration shell model
        autoPtr<baffleType> baffle_;

        //- Dictionary
        dictionary dict_;


public:

    //- Runtime type information
    TypeName("vibrationShell");


    // Constructors

        //- Construct from patch and internal field
        vibrationShellFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        vibrationShellFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- vibrationShellFvPatchScalarField onto a new patch
        vibrationShellFvPatchScalarField
        (
            const vibrationShellFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        vibrationShellFvPatchScalarField
        (
            const vibrationShellFvPatchScalarField&,
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


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
