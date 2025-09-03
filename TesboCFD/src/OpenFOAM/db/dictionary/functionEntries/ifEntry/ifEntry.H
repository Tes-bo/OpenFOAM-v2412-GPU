/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::functionEntries::ifEntry

Description
    Conditional parsing of dictionary entries.

    E.g.
    \verbatim
        U_inlet 15;

        #if #calc "${U_inlet} < 10"
            ..
        #else
            ..
        #endif
    \endverbatim

    Note:
    - only supports single line, '\' is not supported
    - condition should be readable as a \c Switch
      (0,1, true, false, etc.) or a scalar (0.0, ...)

See also
    Foam::functionEntries::ifeqEntry
    Foam::Switch

SourceFiles
    ifEntry.C

\*---------------------------------------------------------------------------*/

#ifndef functionEntries_ifEntry_H
#define functionEntries_ifEntry_H

#include "ifeqEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                          Class ifEntry Declaration
\*---------------------------------------------------------------------------*/

class ifEntry
:
    public ifeqEntry
{
    friend class ifeqEntry;

    // Private Member Functions

        //- Test first token (label, scalar, word) for true/false
        static bool isTrue(ITstream& is);

        //- Execute the functionEntry in a sub-dict context
        static bool execute
        (
            DynamicList<filePos>& stack,
            dictionary& parentDict,
            Istream& is
        );

public:

    //- Runtime type information
    ClassName("if");

    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
