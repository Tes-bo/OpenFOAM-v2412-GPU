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
    Foam::uniformFixedGradientFaPatchField

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
    The uniformGradient entry is a Function1 type.
    The example above gives the usage for supplying a constant value.

    The \c value entry (optional) is used for the initial values.
    Otherwise the \c uniformGradient is used for the evaluation.
    In some cases (eg, coded or expression entries with references to other
    fields) this can be problematic and the \c value entry will be needed.

See also
    Foam::Function1Types
    Foam::fixedGradientFaPatchField

SourceFiles
    uniformFixedGradientFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_uniformFixedGradientFaPatchField_H
#define Foam_uniformFixedGradientFaPatchField_H

#include "fixedGradientFaPatchField.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class uniformFixedGradientFaPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformFixedGradientFaPatchField
:
    public fixedGradientFaPatchField<Type>
{
    // Private Data

        //- Function providing the gradient
        autoPtr<Function1<Type>> refGradFunc_;


public:

    //- Runtime type information
    TypeName("uniformFixedGradient");


    // Constructors

        //- Construct from patch and internal field
        uniformFixedGradientFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch and internal field and patch field
        uniformFixedGradientFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const Field<Type>& fld
        );

        //- Construct from patch, internal field and dictionary
        uniformFixedGradientFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        uniformFixedGradientFaPatchField
        (
            const uniformFixedGradientFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        uniformFixedGradientFaPatchField
        (
            const uniformFixedGradientFaPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        uniformFixedGradientFaPatchField
        (
            const uniformFixedGradientFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return clone
        virtual tmp<faPatchField<Type>> clone() const
        {
            return faPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return faPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write includes "value" entry
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uniformFixedGradientFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
