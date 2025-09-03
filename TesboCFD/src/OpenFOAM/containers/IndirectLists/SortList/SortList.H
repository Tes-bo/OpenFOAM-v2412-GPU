/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::SortList

Description
    An indirect list with addressing based on sorting.
    The list is sorted upon construction or when explicitly requested.

    Uses the std::stable_sort() algorithm.

SourceFiles
    SortListI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_SortList_H
#define Foam_SortList_H

#include "IndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class SortList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class SortList
:
    public IndirectList<T>
{
public:

    // Constructors

        //- Shallow copy values list reference, sort immediately
        inline explicit SortList(const UList<T>& values);

        //- Shallow copy values list reference,
        //- sort with given \b value comparator.
        //  \note The comparator is not stored in the class.
        template<class Compare>
        inline SortList(const UList<T>& values, const Compare& comp);


    // Member Functions

        //- Return the list of sorted indices (updated every sort).
        //  Same as addressing()
        inline const labelUList& indices() const noexcept;

        //- Return the list of indices (updated every sort).
        //  Same as addressing()
        inline labelList& indices() noexcept;

        //- Reverse the indices
        inline void reverse();

        //- Reset list indices to identity
        inline void reset();

        //- Sort the list using the given \b value comparator
        template<class Compare>
        inline void sort(const Compare& comp);

        //- Forward (stable) sort the list.
        //- Functionally identical to sort with std::less\<T\>()
        inline void sort();

        //- Reverse (stable) sort the list.
        //- Functionally identical to sort with std::greater\<T\>()
        inline void reverseSort();

        //- Sort the list, only retaining unique entries
        inline void uniqueSort();


    // Member Operators

        //- Assignment operators
        using IndirectList<T>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SortListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
