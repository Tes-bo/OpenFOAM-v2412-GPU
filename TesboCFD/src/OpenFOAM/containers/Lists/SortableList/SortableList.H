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
    Foam::SortableList

Description
    A list that is sorted upon construction or when explicitly requested
    with the sort() method.

    Uses the std::stable_sort() algorithm.

Note
    In many cases you may wish to reuse list storage.
    The Foam::sortedOrder() function and the Foam::SortList container
    provide two other alternatives.

SourceFiles
    SortableList.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SortableList_H
#define Foam_SortableList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SortableList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class SortableList
:
    public List<T>
{
    // Private Data

        //- Indices from last sort()
        labelList indices_;


public:

    // Constructors

        //- Default construct
        SortableList() noexcept = default;

        //- Construct given size, sort later.
        //  The indices remain empty until the list is sorted
        inline explicit SortableList(const label size);

        //- Construct zero-initialized with given size, sort later.
        //  The indices remain empty until the list is sorted
        inline SortableList(const label size, const Foam::zero);

        //- Construct given size and initial value, sorting later.
        //  The indices remain empty until the list is sorted
        inline SortableList(const label size, const T& val);

        //- Copy construct
        inline SortableList(const SortableList<T>& lst);

        //- Move construct
        inline SortableList(SortableList<T>&& lst);

        //- Copy construct from UList, sorting immediately
        explicit inline SortableList(const UList<T>& values);

        //- Move construct from List, sorting immediately
        inline SortableList(List<T>&& values);

        //- Construct from an initializer list, sorting immediately
        inline SortableList(std::initializer_list<T> values);


    // Member Functions

        //- Return the list of sorted indices. Updated every sort
        const labelList& indices() const noexcept
        {
            return indices_;
        }

        //- Return non-const access to the sorted indices. Updated every sort
        labelList& indices() noexcept
        {
            return indices_;
        }

        //- Clear the list and the indices
        inline void clear();

        //- Clear the indices and return a reference to the underlying List
        inline List<T>& shrink();

        //- Forward (stable) sort the list (if changed after construction).
        //  Resizes the indices as required
        void sort();

        //- Reverse (stable) sort the list
        //  Resizes the indices as required
        void reverseSort();

        //- Forward partial sort the list until the middle point
        void partialSort(label n, label start=0);

        //- Reverse partial sort the list until the middle point
        void partialReverseSort(label n, label start=0);

        //- Swap content with another SortableList in constant time
        inline void swap(SortableList<T>& other);


    // Member Operators

        //- Assignment of all entries to the given value, removing indices.
        inline void operator=(const T& val);

        //- Assignment to UList operator, removing indices. Takes linear time
        inline void operator=(const UList<T>& lst);

        //- Assignment operator. Takes linear time
        inline void operator=(const SortableList<T>& lst);

        //- Move assignment, removing indices. Takes linear time
        inline void operator=(List<T>&& lst);

        //- Move operator (constant time)
        inline void operator=(SortableList<T>&& lst);

        //- Assignment to an initializer list
        void operator=(std::initializer_list<T> lst);
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Does not need std::swap or Foam::Swap() specialization
// since SortableList is MoveConstructible and MoveAssignable


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SortableList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
