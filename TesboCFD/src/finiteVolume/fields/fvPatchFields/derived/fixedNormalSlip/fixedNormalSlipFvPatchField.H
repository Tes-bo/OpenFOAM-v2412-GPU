/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fixedNormalSlipFvPatchField

Group
    grpGenericBoundaryConditions grpWallBoundaryConditions

Description
    This boundary condition sets the patch-normal component to the field (vector
    or tensor) to the patch-normal component of a user specified field.
    The tangential component is treated as slip, i.e. copied from the internal
    field.

    The "value" entry is NO_READ, optional write.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type            fixedNormalSlip;
        fixedValue      uniform (1 0 0);    // example entry for a vector field

        // Optional entries
        writeValue      false;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                      | Type    | Reqd | Deflt
      type         | Type name: fixedNormalSlip       | word    | yes  | -
      fixedValue   | User-defined value the normal component of which  <!--
                   --> the boundary is set to                          <!--
                   -->                | vectorField/tensorField | yes  | -
      writeValue   | Output patch values (eg, ParaView) | bool  | no   | false
    \endtable

    The inherited entries are elaborated in:
      - \link transformFvPatchField.H \endlink
      - \link fvPatchField.H \endlink

See also
    Foam::transformFvPatchField

SourceFiles
    fixedNormalSlipFvPatchField.C
    fixedNormalSlipFvPatchFields.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fixedNormalSlipFvPatchField_H
#define Foam_fixedNormalSlipFvPatchField_H

#include "transformFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fixedNormalSlipFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedNormalSlipFvPatchField
:
    public transformFvPatchField<Type>
{
    //- The parent boundary condition type
    typedef transformFvPatchField<Type> parent_bctype;


    // Private Data

        //- Value the normal component of which the boundary is set to
        Field<Type> fixedValue_;

        //- Flag to output patch values (e.g. for ParaView)
        bool writeValue_;


public:

    //- Runtime type information
    TypeName("fixedNormalSlip");


    // Constructors

        //- Construct from patch and internal field
        fixedNormalSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary.
        //- The "value" entry is NO_READ.
        fixedNormalSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedNormalSlipFvPatchField
        //- onto a new patch
        fixedNormalSlipFvPatchField
        (
            const fixedNormalSlipFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedNormalSlipFvPatchField
        (
            const fixedNormalSlipFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        fixedNormalSlipFvPatchField
        (
            const fixedNormalSlipFvPatchField<Type>&,
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

        //- False: this patch field is not altered by assignment.
        virtual bool assignable() const { return false; }


        // Access

            //- User-defined input field (modifiable field)
            virtual Field<Type>& fixedValue()
            {
                return fixedValue_;
            }

            //- User-defined input field (const field)
            virtual const Field<Type>& fixedValue() const
            {
                return fixedValue_;
            }


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper& m
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>& ptf,
                const labelList& addr
            );


        // Evaluation

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::buffered
            );

            //- Return face-gradient transform diagonal
            virtual tmp<Field<Type>> snGradTransformDiag() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fixedNormalSlipFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
