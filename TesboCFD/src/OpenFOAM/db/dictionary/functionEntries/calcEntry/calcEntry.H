/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::functionEntries::calcEntry

Description
    Uses dynamic compilation to provide calculating functionality
    for entering dictionary entries.

    E.g.

    \verbatim
    a 1.0;
    b 3;
    c \#calc "$a/$b";
    \endverbatim

    Note the explicit trailing 0 ('1.0') to force a to be read (and written)
    as a floating point number.

Note
    Internally this is just a wrapper around codeStream functionality - the
    \#calc string gets used to construct a dictionary for codeStream.

SourceFiles
    calcEntry.C

\*---------------------------------------------------------------------------*/

#ifndef functionEntries_calcEntry_H
#define functionEntries_calcEntry_H

#include "codeStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                         Class calcEntry Declaration
\*---------------------------------------------------------------------------*/

class calcEntry
:
    public codeStream
{
    // Private Member Functions

    //- Evaluate dynamically compiled code, returning result as string
    static string evaluate(const dictionary& parentDict, Istream& is);


public:

    //- Execute in a primitiveEntry context
    static bool execute
    (
        const dictionary& parentDict,
        primitiveEntry& entry,
        Istream& is
    );

    //- Execute in a sub-dict context
    static bool execute(dictionary& parentDict, Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
