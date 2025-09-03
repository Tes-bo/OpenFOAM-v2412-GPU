/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

Namespace
    Foam::HashSetOps

Description
    Various operations for HashSet.

Namespace
    Foam::HashTableOps

Description
    Various operations for HashTable.

SourceFiles
    HashOps.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_HashOps_H
#define Foam_HashOps_H

#include "HashSet.H"
#include "List.H"

namespace Foam
{

// Forward Declarations
class bitSet;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       Namespace HashSetOps Declaration
\*---------------------------------------------------------------------------*/

namespace HashSetOps
{

//- Combine HashSet operation. Equivalent to 'a |= b'
template<class Key=word, class HashType=Foam::Hash<Key>>
struct plusEqOp
{
    typedef HashSet<Key, HashType> value_type;

    void operator()(value_type& a, const value_type& b) const
    {
        a |= b;
    }
};


//- Convert a bitset to a labelHashSet of the indices used.
//
//  \param select the bitset for which an \a on entry corresponds
//       to an index in the output labelHashSet
//
//  \return a labelHashSet of the selected indices
//
//  This is equivalent of the following code, but more efficiently implemented.
//  \code
//      bitSet select = ...;
//      return labelHashSet(select.toc());
//  \endcode
labelHashSet used(const bitSet& select);


//- Convert a list of bools to a labelHashSet of the indices used.
//
//  \param select the boolList for which a \a true entry corresponds
//       to an index in the output labelHashSet
//
//  \return a labelHashSet of the selected indices
labelHashSet used(const UList<bool>& select);


//- Transform the \a on locations to a bitSet.
//  Ignored any negative values (invalid positions in a bitset).
//
//  \param locations the list of positions corresponding to an \a on bit.
//
//  \return a bitset
//
//  \see BitSetOps::create for other possibilities
bitSet bitset(const labelHashSet& locations);


//- Transform the \a on locations to a boolList, with \a true for each
//- non-negative location and \a false for all others.
//
//  \param locations the list of positions corresponding to an \a on bit.
//
//  \return a boolList
//
//  \note The operation necessarily discards any negative values since these
//     are invalid positions in a boolList.
List<bool> bools(const labelHashSet& locations);

} // End namespace HashSetOps


/*---------------------------------------------------------------------------*\
                      Namespace HashTableOps Declaration
\*---------------------------------------------------------------------------*/

namespace HashTableOps
{

//- Combine HashTable operation. Equivalent to 'a += b'
template<class T, class Key=word, class HashType=Foam::Hash<Key>>
struct plusEqOp
{
    typedef HashTable<T, Key, HashType> value_type;

    void operator()(value_type& a, const value_type& b) const
    {
        a += b;
    }
};


//- List of values from HashTable, optionally sorted.
template<class T, class Key, class Hash>
List<T> values
(
    const HashTable<T, Key, Hash>& tbl,
    const bool doSort=false
)
{
    List<T> output(tbl.size());

    label i=0;

    forAllConstIters(tbl, iter)
    {
        output[i] = iter.val();
        ++i;
    }

    if (doSort)
    {
        Foam::sort(output);
    }

    return output;
}

} // End namespace HashTableOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
