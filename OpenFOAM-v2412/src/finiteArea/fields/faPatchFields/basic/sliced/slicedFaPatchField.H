/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::slicedFaPatchField

Group
    grpGenericBoundaryConditions

Description
    Specialization of faPatchField which creates the underlying
    faPatchField as a slice of the given complete field.

    The destructor is wrapped to avoid deallocation of the storage of the
    complete fields when this is destroyed.

    Should only used as a template argument for SlicedGeometricField.

See also
    Foam::faPatchField

SourceFiles
    slicedFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_slicedFaPatchField_H
#define Foam_slicedFaPatchField_H

#include "faPatchField.H"
#include "processorFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class slicedFaPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class slicedFaPatchField
:
    public faPatchField<Type>
{
public:

    //- Type for the mesh processor patch
    typedef processorFaPatch processorPatchType;


    //- Runtime type information
    TypeName("sliced");


    // Constructors

        //- Construct from patch, internal field and field to slice
        slicedFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const Field<Type>& completeOrBoundaryField,
            const bool isBoundaryOnly = false
        );

        //- Construct from patch and internal field. Assign value later.
        slicedFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary.
        //- Not implemented
        slicedFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given patch field onto a new patch.
        //- Not implemented
        slicedFaPatchField
        (
            const slicedFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        slicedFaPatchField(const slicedFaPatchField<Type>&);

        //- Construct as copy setting internal field reference
        slicedFaPatchField
        (
            const slicedFaPatchField<Type>&,
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


    //- Destructor, resets pointers to avoid deletion of underlying field
    virtual ~slicedFaPatchField();


    // Member Functions

        //- True: this patch field fixes a value.
        virtual bool fixesValue() const { return true; }

        //- False: this patch field is not altered by assignment.
        virtual bool assignable() const { return false; }


        // Evaluation functions

            //- Return patch-normal gradient
            virtual tmp<Field<Type>> snGrad() const;

            //- Return internal field next to patch as patch field
            virtual tmp<Field<Type>> patchInternalField() const;

            //- Return internal field next to patch as patch field
            virtual void patchInternalField(Field<Type>&) const;

            //- Return neighbour coupled given internal cell data
            virtual tmp<Field<Type>> patchNeighbourField
            (
                const Field<Type>& iField
            ) const;

            //- Return patchField of the values on the patch or on the
            //  opposite patch
            virtual tmp<Field<Type>> patchNeighbourField() const;

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::buffered
            )
            {}

            //- Evaluate the patch field, sets updated() to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::buffered
            )
            {}

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const faPatchField<Type>&) {}
        virtual void operator+=(const faPatchField<Type>&) {}
        virtual void operator-=(const faPatchField<Type>&) {}
        virtual void operator*=(const faPatchField<scalar>&) {}
        virtual void operator/=(const faPatchField<scalar>&) {}

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
    #include "slicedFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
