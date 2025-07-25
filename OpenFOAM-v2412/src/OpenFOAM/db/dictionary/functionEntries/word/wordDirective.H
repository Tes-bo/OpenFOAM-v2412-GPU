/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::functionEntries::wordDirective

Description
    Converts/expands string content to a word.
    This can be useful for composition of names.

    For example,
    \verbatim
    #word "some_prefix_solverInfo_${application}_${{10 * 5}}"
    {
        type    solverInfo;
        libs    (utilityFunctionObjects);
        fields  (".*");

        returnName  #word solver$type;
    }
    \endverbatim

SourceFiles
    wordDirective.C

\*---------------------------------------------------------------------------*/

#ifndef functionEntries_wordDirective_H
#define functionEntries_wordDirective_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class wordDirective Declaration
\*---------------------------------------------------------------------------*/

class wordDirective
:
    public functionEntry
{
    // Private Member Functions

    //- Evaluate
    static token evaluate(const dictionary& parentDict, Istream& is);


public:

    //- Execute in a primitiveEntry context
    static bool execute
    (
        const dictionary& parentDict,
        primitiveEntry& entry,
        Istream& is
    );

    //- Execute in a sub-dict context.
    static bool execute(dictionary& parentDict, Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
