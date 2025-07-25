/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::IjkField

Description
    Generic templated field type with i-j-k addressing.

SourceFiles
    IjkFieldI.H
    IjkField.C

\*---------------------------------------------------------------------------*/

#ifndef IjkField_H
#define IjkField_H

#include "ijkAddressing.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IjkField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IjkField
:
    public Field<Type>
{
    // Private Data

        //- The i-j-k addressing
        ijkAddressing ijk_;


public:

    // Constructors

        //- Construct zero-sized
        inline IjkField();

        //- Copy construct
        inline IjkField(const IjkField<Type>& field);

        //- Move construct
        inline IjkField(IjkField<Type>&& field);

        //- Construct with sizing information, leaving values uninitialized
        inline explicit IjkField(const labelVector& ijk);

        //- Construct with sizing information and initial value
        inline IjkField(const labelVector& ijk, const Type& val);

        //- Construct with sizing information and initial values of zero
        inline IjkField(const labelVector& ijk, const Foam::zero);

        //- Copy construct from components
        inline IjkField(const labelVector& ijk, const UList<Type>& list);

        //- Move construct from components
        inline IjkField(const labelVector& ijk, Field<Type>&& field);


    // Member Functions

    // Access

        //- Return i,j,k addressing
        inline const ijkAddressing& ijk() const;

        //- Return i,j,k addressing for modification
        inline ijkAddressing& ijk();

        //- Return i,j,k addressing sizes
        inline const labelVector& sizes() const;

        //- Return i,j,k addressing sizes for modification
        inline labelVector& sizes();

        //- The field size
        using Field<Type>::size;

        //- The addressing dimension in the given direction
        inline const label& size(const vector::components cmpt) const;


    // Edit

        //- Clear dimensions and field
        inline void clear();

        //- Change dimensions. Fill new values with Zero
        void resize(const labelVector& newSizes);

        //- Change dimensions
        void resize(const labelVector& newSizes, const Type& val);


    // Access Operators

        //- Field access at given i-j-k position
        inline const Type& operator()
        (
            const label i,
            const label j,
            const label k
        ) const;

        //- Field access at given i-j-k position
        inline Type& operator()
        (
            const label i,
            const label j,
            const label k
        );

        //- Field access at given i-j-k position
        inline const Type& operator()(const labelVector& ijk) const;

        //- Field access at given i-j-k position
        inline Type& operator()(const labelVector& ijk);


    // Member Operators

        //- Copy assignment
        void operator=(const IjkField<Type>& rhs);
        void operator=(const tmp<IjkField<Type>>& rhs);

        //- Move assignment
        inline void operator=(IjkField<Type>&& rhs);

        //- Value assignment
        inline void operator=(const Type& val);
        inline void operator=(const Foam::zero);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IjkFieldI.H"

#ifdef NoRepository
    #include "IjkField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
