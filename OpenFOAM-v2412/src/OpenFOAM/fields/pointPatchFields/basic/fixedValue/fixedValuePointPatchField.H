/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::fixedValuePointPatchField

Description
    A FixedValue boundary condition for pointField.

    The "value" entry is normally MUST_READ.

SourceFiles
    fixedValuePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fixedValuePointPatchField_H
#define Foam_fixedValuePointPatchField_H

#include "valuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fixedValuePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedValuePointPatchField
:
    public valuePointPatchField<Type>
{
public:

    //- Runtime type information
    TypeName("fixedValue");


    // Constructors

        //- Construct from patch and internal field
        fixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and value
        fixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const Type& value
        );

        //- Construct from patch, internal field and dictionary
        fixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&,
            IOobjectOption::readOption requireValue = IOobjectOption::MUST_READ
        );

        //- Compatibility. Prefer with readOption
        fixedValuePointPatchField
        (
            const pointPatch& p,
            const DimensionedField<Type, pointMesh>& iF,
            const dictionary& dict,
            const bool needValue
        )
        :
            fixedValuePointPatchField
            (
                p, iF, dict,
                (needValue? IOobjectOption::MUST_READ : IOobjectOption::NO_READ)
            )
        {}

        //- Construct by mapping given patch field onto a new patch
        fixedValuePointPatchField
        (
            const fixedValuePointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        fixedValuePointPatchField
        (
            const fixedValuePointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const
        {
            return pointPatchField<Type>::Clone(*this);
        }

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return pointPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        //- True: this patch field fixes a value.
        virtual bool fixesValue() const { return true; }


    // Member Operators

        // Disable assignment operators

            virtual void operator=(const Field<Type>&){}
            virtual void operator=(const Type&){}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fixedValuePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
