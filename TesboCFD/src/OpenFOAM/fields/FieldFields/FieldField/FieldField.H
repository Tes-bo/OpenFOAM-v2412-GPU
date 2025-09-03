/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
    Foam::FieldField

Description
    A field of fields is a PtrList of fields with reference counting.

SourceFiles
    FieldField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_FieldField_H
#define Foam_FieldField_H

#include "tmp.H"
#include "Pair.H"
#include "PtrList.H"
#include "scalar.H"
#include "direction.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<template<class> class Field, class Type>
class FieldField;

template<template<class> class Field, class Type>
Ostream& operator<<
(
    Ostream&,
    const FieldField<Field, Type>&
);

template<template<class> class Field, class Type>
Ostream& operator<<
(
    Ostream&,
    const tmp<FieldField<Field, Type>>&
);


/*---------------------------------------------------------------------------*\
                           Class FieldField Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class Field, class Type>
class FieldField
:
    public refCount,
    public PtrList<Field<Type>>
{
public:

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;


    // Constructors

        //- Construct null
        //  Used for temporary fields which are initialised after construction
        constexpr FieldField() noexcept;

        //- Construct given size
        //  Used for temporary fields which are initialised after construction
        explicit FieldField(const label size);

        //- Clone construct with new type
        FieldField(const word& type, const FieldField<Field, Type>& ff);

        //- Copy construct, cloning each element
        FieldField(const FieldField<Field, Type>& ff);

        //- Move construct
        FieldField(FieldField<Field, Type>&& ff);

        //- Construct as copy or re-use as specified.
        FieldField(FieldField<Field, Type>& ff, bool reuse);

        //- Copy construct from PtrList
        FieldField(const PtrList<Field<Type>>& list);

        //- Move construct from PtrList
        FieldField(PtrList<Field<Type>>&& list);

        //- Move/copy construct from tmp<FieldField>
        FieldField(const tmp<FieldField<Field, Type>>& tf);

        //- Construct from Istream
        FieldField(Istream& is);

        //- Clone
        tmp<FieldField<Field, Type>> clone() const;

        //- Return a pointer to a new calculatedFvPatchFieldField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<FieldField<Field, Type>> NewCalculatedType
        (
            const FieldField<Field, Type2>& ff
        );


    // Member Functions

        //- Negate this field. See notes in Field
        void negate();

        //- Normalise this field. See notes in Field
        void normalise();

        //- Return a component field of the field
        tmp<FieldField<Field, cmptType>> component(const direction) const;

        //- Replace a component field of the field
        void replace(const direction, const FieldField<Field, cmptType>&);

        //- Replace a component field of the field
        void replace(const direction, const cmptType&);

        //- Impose lower (floor) clamp on the field values (in-place)
        void clamp_min(const Type& lower);

        //- Impose upper (ceiling) clamp on the field values (in-place)
        void clamp_max(const Type& upper);

        //- Clamp field values (in-place) to the specified range.
        //  Does not check if range is valid or not.
        void clamp_range(const Type& lower, const Type& upper);

        //- Clamp field values (in-place) to the specified range.
        //  Does not check if range is valid or not.
        void clamp_range(const MinMax<Type>& range);

        //- Return the field transpose (only defined for second rank tensors)
        tmp<FieldField<Field, Type>> T() const;


    // Member Operators

        //- Const or non-const access to a field
        using PtrList<Field<Type>>::operator[];

        //- Const access to a single field element via (fieldi, elemi)
        inline const Type& operator[](const labelPair& index) const;

        //- Non-const access to a single field element via (fieldi, elemi)
        inline Type& operator[](const labelPair& index);

        //- Copy assignment
        void operator=(const FieldField<Field, Type>&);

        //- Move assignment
        void operator=(FieldField<Field, Type>&&);

        //- Move or clone assignment
        void operator=(const tmp<FieldField<Field, Type>>&);

        //- Assign uniform value
        void operator=(const Type& val);

        //- Assign uniform zero
        void operator=(const Foam::zero);

        void operator+=(const FieldField<Field, Type>&);
        void operator+=(const tmp<FieldField<Field, Type>>&);

        void operator-=(const FieldField<Field, Type>&);
        void operator-=(const tmp<FieldField<Field, Type>>&);

        void operator*=(const FieldField<Field, scalar>&);
        void operator*=(const tmp<FieldField<Field, scalar>>&);

        void operator/=(const FieldField<Field, scalar>&);
        void operator/=(const tmp<FieldField<Field, scalar>>&);

        void operator+=(const Type&);
        void operator-=(const Type&);

        void operator*=(const scalar&);
        void operator/=(const scalar&);


    // IOstream operators

        friend Ostream& operator<< <Field, Type>
        (
            Ostream&,
            const FieldField<Field, Type>&
        );

        friend Ostream& operator<< <Field, Type>
        (
            Ostream&,
            const tmp<FieldField<Field, Type>>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldFieldFunctions.H"

#ifdef NoRepository
    #include "FieldField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
