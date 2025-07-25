/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::IndirectSubList

Description
    Indirect access to a sub-section of a list.

    In many cases, using a Foam::SubList provides a simpler and more
    efficient means of accessing a sub-list.
    There are, however, some advantages of a IndirectSubList:
    - allows adjustment of its addressing range after construct
    - can recover the original, underlying list at any time

SeeAlso
    Foam::SubList

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef IndirectSubList_H
#define IndirectSubList_H

#include "IndirectListAddressing.H"
#include "IndirectListBase.H"
#include "labelRange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class IndirectSubList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IndirectSubList
:
    private IndirectListAddressing<labelRange>,
    public IndirectListBase<T, labelRange>
{
public:

    // Constructors

        //- Construct from UList, the entire size
        explicit IndirectSubList(const UList<T>& values)
        :
            IndirectListAddressing<labelRange>(labelRange(0, values.size())),
            IndirectListBase<T, labelRange>
            (
                values,
                IndirectListAddressing<labelRange>::addressing()
            )
        {}

        //- Construct from values list and range
        IndirectSubList
        (
            const UList<T>& values,
            const labelRange& addr
        )
        :
            IndirectListAddressing<labelRange>(addr),
            IndirectListBase<T, labelRange>
            (
                values,
                IndirectListAddressing<labelRange>::addressing()
            )
        {}


    // Member Functions

        //- The list addressing
        using IndirectListAddressing::addressing;


    // Member Operators

        //- Use standard assignment operations
        using IndirectListBase<T, labelRange>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
