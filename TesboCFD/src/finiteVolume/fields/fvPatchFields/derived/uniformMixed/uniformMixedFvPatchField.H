/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::uniformMixedFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides 'mixed' type boundary condition
    that mix a \em uniform fixed value and a \em uniform patch-normal
    gradient condition. The term "uniform" is a legacy name since the
    prescribed values were previously spatially uniform across that patch.

    In the meantime, a PatchFunction1 is used, which can have both
    spatial and temporal dependencies.

Usage
    \table
        Property     | Description                      | Required | Default
        uniformValue | uniform value                    | partly | 0
        uniformGradient | uniform gradient              | partly | 0
        uniformValueFraction | uniform valueFraction    | partly | depends
        value        | initial field value              | optional |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformMixed;
        uniformValue    constant 0.2;
        uniformGradient constant 0.2;
        uniformValueFraction
        {
            type  sine;
            ...
        }
    }
    \endverbatim

Note
    This boundary condition allows \em lazier definitions so that either
    or both: \c uniformValue and \c uniformGradient must be defined.
    If only of these entries is defined, the value fraction is automatically
    treated appropriately (ie, 0 with \c uniformGradient and 1 with
    uniformValue).
    If both \c uniformValue and \c uniformGradient are defined,
    the \c uniformValueFraction must also be defined.

    The \c value entry (optional) is used for the initial values.
    Otherwise the function(s) are used for the evaluation.
    In some cases (eg, coded or expression entries with references to other
    fields) this can be problematic and the \c value entry will be needed.

See also
    Foam::Function1Types
    Foam::mixedFvPatchField

SourceFiles
    uniformMixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_uniformMixedFvPatchField_H
#define Foam_uniformMixedFvPatchField_H

#include "mixedFvPatchField.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class uniformMixedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformMixedFvPatchField
:
    public mixedFvPatchField<Type>
{
    // Private Data

        //- Function providing the value
        autoPtr<PatchFunction1<Type>> refValueFunc_;

        //- Function providing the gradient
        autoPtr<PatchFunction1<Type>> refGradFunc_;

        //- Function providing the value-fraction
        autoPtr<PatchFunction1<scalar>> valueFractionFunc_;


public:

    //- Runtime type information
    TypeName("uniformMixed");


    // Constructors

        //- Construct from patch and internal field
        uniformMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        uniformMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>& fld
        );

        //- Construct from patch, internal field and dictionary
        uniformMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        uniformMixedFvPatchField
        (
            const uniformMixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformMixedFvPatchField
        (
            const uniformMixedFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        uniformMixedFvPatchField
        (
            const uniformMixedFvPatchField<Type>&,
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


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uniformMixedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
