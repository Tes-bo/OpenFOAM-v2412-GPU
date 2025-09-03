/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::instant

Description
    An instant of time. Contains the time value and name.
    Uses Foam::Time when formatting the name.

SourceFiles
    instant.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_instant_H
#define Foam_instant_H

#include "Instant.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class Pair;
template<class T> class UList;

/*---------------------------------------------------------------------------*\
                           Class instant Declaration
\*---------------------------------------------------------------------------*/

class instant
:
    public Instant<word>
{
public:

    // Static Data Members

        //- The type name (eg, for pTraits)
        static const char* const typeName;


    // Generated Methods

        //- Default construct, with value = 0 and empty name
        instant() = default;

        //- Copy construct
        instant(const instant&) = default;

        //- Move construct
        instant(instant&&) = default;

        //- Copy assignment
        instant& operator=(const instant&) = default;

        //- Move assignment
        instant& operator=(instant&&) = default;


    // Constructors

        //- Inherit all constructors (eg, copy/move construct from components)
        using Instant<word>::Instant;

        //- Construct from timeValue, auto generating the name
        explicit instant(scalar timeValue);

        //- Copy construct from timeName, parsing timeName for a value
        explicit instant(const word& timeName);

        //- Move construct from timeName, parsing timeName for a value
        explicit instant(word&& timeName);


    // Helper Functions (for searching)

        //- Find and return index of given start time (linear search)
        static label findStart
        (
            const UList<instant>& times,
            const scalar timeVal
        );

        //- Find lower/upper indices for given time value in list of instances
        //- (linear search) continuing \em after the given start index.
        //  \returns the range indices or (-1,-1) if unsuccessful.
        static Pair<label> findRange
        (
            const UList<instant>& times,
            const scalar timeVal,
            const label start = -1
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
