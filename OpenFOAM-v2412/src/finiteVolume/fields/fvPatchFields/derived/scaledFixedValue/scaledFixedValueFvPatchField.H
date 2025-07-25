/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::scaledFixedValueFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This condition applies a scalar multiplier to the value of another
    boundary condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        scale        | Time varying scale      | yes         |
        refValue     | patchField providing the raw patch value | yes |
    \endtable

    Example of the boundary condition specification to scale a reference
    velocity of (15 0 0)  supplied as a fixedValue by a table of values
    that ramps the scale from 0 to 1 over 1 second:
    \verbatim
    <patchName>
    {
        type            scaledFixedValue;

        scale table
        (
            (    0   0)
            (  1.0 1.0)
            (100.0 1.0)
        );

        refValue
        {
            type            fixedValue;
            value           uniform (15 0 0);
        }
    }
    \endverbatim

SourceFiles
    scaledFixedValueFvPatchField.C

SeeAlso
    Foam::PatchFunction1

\*---------------------------------------------------------------------------*/

#ifndef scaledFixedValueFvPatchField_H
#define scaledFixedValueFvPatchField_H

#include "fixedValueFvPatchField.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class scaledFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class scaledFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
protected:

    // Protected data

        //- Scalar scale factor
        autoPtr<PatchFunction1<scalar>> scalePtr_;

        //- Condition to supply the reference value
        tmp<fvPatchField<Type>> refValuePtr_;


public:

    //- Runtime type information
    TypeName("scaledFixedValue");


    // Constructors

        //- Construct from patch and internal field
        scaledFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        scaledFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a scaledFixedValueFvPatchField onto
        //- a new patch
        scaledFixedValueFvPatchField
        (
            const scaledFixedValueFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        scaledFixedValueFvPatchField(const scaledFixedValueFvPatchField&);

        //- Construct as copy setting internal field reference
        scaledFixedValueFvPatchField
        (
            const scaledFixedValueFvPatchField&,
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

        // Access

            //- Return the reference value condition
            const fvPatchField<Type>& refValue() const
            {
                return refValuePtr_();
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper& m);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>& ptf,
                const labelList& addr
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        // Force an assignment, applying scaling to refValue

        virtual void operator==(const fvPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);

        // Prevent automatic comparison rewriting (c++20)
        bool operator!=(const fvPatchField<Type>&) const = delete;
        bool operator!=(const Field<Type>&) const = delete;
        bool operator!=(const Type&) const = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "scaledFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
