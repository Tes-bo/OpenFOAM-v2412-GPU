/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::SliceList

Description
    A List with indirect slice addressing.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef SliceList_H
#define SliceList_H

#include "IndirectListAddressing.H"
#include "IndirectListBase.H"
#include "sliceRange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class SliceList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class SliceList
:
    private IndirectListAddressing<sliceRange>,
    public IndirectListBase<T, sliceRange>
{
public:

    // Constructors

        //- Copy construct from values list and slicing
        SliceList(const UList<T>& values, const sliceRange& addr)
        :
            IndirectListAddressing<sliceRange>(addr),
            IndirectListBase<T, sliceRange>
            (
                values,
                IndirectListAddressing<sliceRange>::addressing()
            )
        {}


    // Member Functions

        //- The list addressing
        using IndirectListAddressing::addressing;


    // Member Operators

        //- Use standard assignment operations
        using IndirectListBase<T, sliceRange>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
