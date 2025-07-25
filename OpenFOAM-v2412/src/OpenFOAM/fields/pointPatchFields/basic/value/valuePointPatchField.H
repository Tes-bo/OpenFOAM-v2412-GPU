/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::valuePointPatchField

Description
    Foam::valuePointPatchField

SourceFiles
    valuePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_valuePointPatchField_H
#define Foam_valuePointPatchField_H

#include "pointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class valuePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class valuePointPatchField
:
    public pointPatchField<Type>,
    public Field<Type>
{
protected:

    // Protected Member Functions

        //- Read the "value" entry into \c *this.
        //  The reading can be optional (default), mandatory etc.
        //  \returns True on success
        bool readValueEntry
        (
            const dictionary& dict,
            IOobjectOption::readOption readOpt = IOobjectOption::LAZY_READ
        );

        //- Write \c *this field as a "value" entry
        void writeValueEntry(Ostream& os) const
        {
            Field<Type>::writeEntry("value", os);
        }

        //- Assign the patch field from the internal field
        void extrapolateInternal();


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("value");


    // Generated Methods

        //- Copy construct
        valuePointPatchField(const valuePointPatchField&) = default;


    // Constructors

        //- Construct from patch and internal field
        valuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and value
        valuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const Type& value
        );

        //- Construct from patch, internal field and dictionary
        valuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary& dict,
            //! The "value" entry (default: mandatory)
            IOobjectOption::readOption requireValue = IOobjectOption::MUST_READ
        );

        //- Construct, forwarding to readOption variant
        valuePointPatchField
        (
            const pointPatch& p,
            const DimensionedField<Type, pointMesh>& iF,
            const dictionary& dict,
            const bool needValue
        )
        :
            valuePointPatchField
            (
                p, iF, dict,
                (needValue? IOobjectOption::MUST_READ : IOobjectOption::NO_READ)
            )
        {}

        //- Construct by mapping given patch field onto a new patch
        valuePointPatchField
        (
            const valuePointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        valuePointPatchField
        (
            const valuePointPatchField<Type>&,
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

        // Access

            //- Return size
            label size() const
            {
                return Field<Type>::size();
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given PointPatchField onto
            //  this PointPatchField
            virtual void rmap
            (
                const pointPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::buffered
            );


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        // Assignment operators

        virtual void operator=(const valuePointPatchField<Type>&);
        virtual void operator=(const pointPatchField<Type>&);
        virtual void operator=(const Field<Type>&);
        virtual void operator=(const Type&);

        // Force an assignment irrespective of form of patch

        virtual void operator==(const valuePointPatchField<Type>&);
        virtual void operator==(const pointPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);

        // Prevent automatic comparison rewriting (c++20)
        bool operator!=(const valuePointPatchField<Type>&) const = delete;
        bool operator!=(const pointPatchField<Type>&) const = delete;
        bool operator!=(const Field<Type>&) const = delete;
        bool operator!=(const Type&) const = delete;
};


// This function is added to override the hack in pointPatchField.H
// which enables simple backward compatibility with versions using
// referenceLevel in GeometricField
template<class Type>
tmp<Field<Type>> operator+
(
    const valuePointPatchField<Type>& vsppf,
    const Type& t
)
{
    return static_cast<const Field<Type>&>(vsppf) + t;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "valuePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
