/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::IndirectList

Description
    A List with indirect addressing.

See also
    Foam::UIndirectList for a version without addressing allocation.

Class
    Foam::UIndirectList

Description
    A List with indirect addressing.
    Like IndirectList but does not store addressing

    Note the const_cast of the list values. This is so we can use it both
    on const and non-const lists. Alternative would be to have a const_
    variant etc.

SourceFiles
    IndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_IndirectList_H
#define Foam_IndirectList_H

#include "List.H"
#include "IndirectListAddressing.H"
#include "IndirectListBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class IndirectList;
template<class T> class UIndirectList;

// Common list types
typedef UIndirectList<bool> boolUIndList;       //!< UIndirectList of bools
typedef UIndirectList<label> labelUIndList;     //!< UIndirectList of labels


/*---------------------------------------------------------------------------*\
                        Class UIndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UIndirectList
:
    public IndirectListBase<T, labelUList>
{
public:

    // Constructors

        //- Shallow copy values and addressing
        UIndirectList(const UList<T>& values, const labelUList& addr)
        :
            IndirectListBase<T, labelUList>(values, addr)
        {}

        //- Copy construct (shallow copy values and addressing)
        UIndirectList(const UIndirectList<T>& list)
        :
            UIndirectList<T>(list.values(), list.addressing())
        {}


    // Member Operators

        //- Use standard assignment operations
        using IndirectListBase<T, labelUList>::operator=;

        //- Deep copy values, Fatal if list sizes are not identical
        void operator=(const UIndirectList<T>& rhs)
        {
            this->copyList(rhs);
        }
};


/*---------------------------------------------------------------------------*\
                        Class IndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IndirectList
:
    private IndirectListAddressing<labelList>,
    public UIndirectList<T>
{
public:

    // Constructors

        //- Copy construct addressing, shallow copy values reference
        inline IndirectList(const UList<T>& values, const labelUList& addr);

        //- Move construct addressing, shallow copy values reference
        inline IndirectList(const UList<T>& values, labelList&& addr);

        //- Zero-sized addressing, shallow copy values reference
        inline IndirectList(const UList<T>& values, const Foam::zero);

        //- Copy construct addressing, shallow copy values reference
        inline IndirectList(const IndirectList<T>& list);

        //- Move construct addressing, shallow copy values reference
        inline IndirectList(IndirectList<T>&& list);

        //- Copy construct addressing, shallow copy values reference
        inline explicit IndirectList(const UIndirectList<T>& list);


    // Static Functions

        //- Return an IndirectList comprising entries with \em positions
        //- that satisfy the condition predicate.
        //
        //  The condition predicate can be considered a type of \em mask
        //  for any given position.
        //  A boolList, bitSet, labelRange or labelHashSet all satisfy
        //  the requirements for use as position condition predicates.
        //
        //  \param values  The source list values
        //  \param select  Accept/reject predicate based on \em position.
        //  \param invert  Invert (negate) the selection logic
        template<class UnaryCondition>
        static IndirectList<T> subset
        (
            const UList<T>& values,
            const UnaryCondition& select,
            const bool invert = false
        );

        //- Return an IndirectList comprising entries with \em values
        //- that satisfy the predicate.
        //
        //  \param values  The source list values
        //  \param pred    Predicate used to test the \em value
        //  \param invert  Invert (negate) the selection logic
        template<class UnaryPredicate>
        static IndirectList<T> subset_if
        (
            const UList<T>& values,
            const UnaryPredicate& pred,
            const bool invert = false
        );

        //- Return an IndirectList with duplicate entries filtered out.
        //  Preserves the original input order, unless sorted = true
        //
        //  \param values  The source list values
        //  \param sorted  Retain sorted order instead of original order
        static IndirectList<T> uniq
        (
            const UList<T>& values,
            const bool sorted = false
        );


    // Member Functions

        //- The list addressing
        using IndirectListAddressing::addressing;


    // Member Operators

        //- Assignment operator
        using UIndirectList<T>::operator=;

        //- Deep copy values, Fatal if list sizes are not identical
        void operator=(const IndirectList<T>& rhs)
        {
            this->copyList(rhs);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IndirectListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
