/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::fixedJumpFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition, using the \c cyclic
    condition as a base.

    The jump is specified as a fixed value field, applied as an offset to the
    'owner' patch.

Usage
    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type should be \c cyclic| yes |
        jump         | current jump value      | yes         |
        relax        | under-relaxation factor | no          |
        minJump      | Minimum jump value      | no          |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedJump;
        patchType       cyclic;
        jump            uniform 10;
    }
    \endverbatim

    The above example shows the use of a fixed jump of '10'.

Note
     The underlying \c patchType should be set to \c cyclic

See also
    Foam::jumpCyclicFvPatchField

SourceFiles
    fixedJumpFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedJumpFvPatchField_H
#define fixedJumpFvPatchField_H

#include "jumpCyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fixedJumpFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedJumpFvPatchField
:
    public jumpCyclicFvPatchField<Type>
{

    // Private data

        //- "jump" field
        Field<Type> jump_;

        //- "jump" field at old time level
        Field<Type> jump0_;

        //- Minimum allowable jump value
        Type minJump_;


        //- Under-relaxation factor
        scalar relaxFactor_;

        //- Time index
        label timeIndex_;


public:

    //- Runtime type information
    TypeName("fixedJump");

    // Constructors

        //- Construct from patch and internal field
        fixedJumpFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedJumpFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            const bool needValue = true
        );

        //- Construct by mapping given fixedJumpFvPatchField onto a
        //  new patch
        fixedJumpFvPatchField
        (
            const fixedJumpFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedJumpFvPatchField
        (
            const fixedJumpFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        fixedJumpFvPatchField
        (
            const fixedJumpFvPatchField<Type>&,
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

            //- Set the jump field
            virtual void setJump(const Field<Type>& jump);

            //- Set the jump field (uniform value)
            virtual void setJump(const Type& jump);

            //- Return the "jump" across the patch
            virtual tmp<Field<Type>> jump() const;

            //- Return the old time "jump" across the patch
            virtual tmp<Field<Type>> jump0() const;

            //- Return the under-relaxation factor
            virtual scalar relaxFactor() const;

            //- Return the relaxed "jump" across the patch
            virtual void relax();


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fixedJumpFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
