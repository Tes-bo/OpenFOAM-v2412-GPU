/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::uniformFixedGradientFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a uniform fixed gradient condition.

Usage
    \table
        Property     | Description                  | Required | Default
        uniformGradient | uniform gradient          | yes |
        value        | initial field value          | optional |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformFixedGradient;
        uniformGradient constant 0.2;
    }
    \endverbatim

Note
    The uniformGradient entry is a PatchFunction1 type,
    able to describe time and spatial varying functions.
    The example above gives the usage for supplying a constant value.

    The \c value entry (optional) is used for the initial values.
    Otherwise the \c uniformGradient is used for the evaluation.
    In some cases (eg, coded or expression entries with references to other
    fields) this can be problematic and the \c value entry will be needed.

See also
    Foam::Function1Types
    Foam::fixedGradientFvPatchField

SourceFiles
    uniformFixedGradientFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_uniformFixedGradientFvPatchField_H
#define Foam_uniformFixedGradientFvPatchField_H

#include "fixedGradientFvPatchFields.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class uniformFixedGradientFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformFixedGradientFvPatchField
:
    public fixedGradientFvPatchField<Type>
{
    // Private Data

        //- Function providing the gradient
        autoPtr<PatchFunction1<Type>> refGradFunc_;


public:

    //- Runtime type information
    TypeName("uniformFixedGradient");


    // Constructors

        //- Construct from patch and internal field
        uniformFixedGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        uniformFixedGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>& fld
        );

        //- Construct from patch, internal field and dictionary
        uniformFixedGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        uniformFixedGradientFvPatchField
        (
            const uniformFixedGradientFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformFixedGradientFvPatchField
        (
            const uniformFixedGradientFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        uniformFixedGradientFvPatchField
        (
            const uniformFixedGradientFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uniformFixedGradientFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
