/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::functionEntries::removeEntry

Description
    Remove a dictionary entry.

    The \c \#remove directive takes a list or a single keyType.
    For example,
    \verbatim
        #remove entry0
        #remove ( entry1 entry2 entry3 otherEntry )
        #remove "entry[1-3]"
        #remove ( "entry[1-3]" otherEntry )
        #remove ^dict1.subdict2.entry2
        #remove "/dict1/subdict2/entry1"
    \endverbatim

Note
    Unless otherwise scoped, the removal occurs in the current context.
    To remove from other scopes, a dot-scoped or slash-scoped syntax is
    required. The slash-scoped syntax must be quoted to ensure that it
    is properly parsed.

SourceFiles
    removeEntry.C

\*---------------------------------------------------------------------------*/

#ifndef functionEntries_removeEntry_H
#define functionEntries_removeEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                         Class removeEntry Declaration
\*---------------------------------------------------------------------------*/

class removeEntry
:
    public functionEntry
{
public:

    //- Remove single or multiple entries. Local or scoped entries.
    static bool execute(dictionary& parentDict, Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
