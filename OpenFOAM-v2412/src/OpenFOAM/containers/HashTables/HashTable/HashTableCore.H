/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::HashTableCore

Description
    Template invariant parts of hash table implementation.

SourceFiles
    HashTableCoreI.H
    HashTableCore.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_HashTableCore_H
#define Foam_HashTableCore_H

#include "label.H"
#include "uLabel.H"
#include "className.H"
#include "nullObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class HashTableCore Declaration
\*---------------------------------------------------------------------------*/

//- Bits that are independent of HashTable template parameters.
struct HashTableCore
{
    //- Maximum allowable internal table size. Approximately labelMax/4
    static const label maxTableSize;

    //- Return a canonical (power-of-two) of the requested size.
    static label canonicalSize(const label requested_size);

    //- Declare type-name (with debug switch)
    ClassName("HashTable");

    //- Default construct
    HashTableCore() noexcept = default;

    static_assert
    (
        sizeof(NullObject) >= sizeof(void*),
        "NullObject is too small to reinterpret_cast as HashTable::iterator"
    );


    //- Factory class for creating a begin/end pair for any const iterator.
    template<class IteratorType, class TableType>
    class const_iterator_pair
    {
        const label size_;
        IteratorType iter_;

    public:

        //- Default construct an empty pair
        inline const_iterator_pair();

        //- Construct begin/end pair for table
        inline const_iterator_pair(const TableType& tbl);

        bool empty() const noexcept { return !size_; }
        label size() const noexcept { return size_; }

        inline IteratorType begin() const;
        inline IteratorType cbegin() const;

        inline IteratorType end() const;
        inline IteratorType cend() const;
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HashTableCoreI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
