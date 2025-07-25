/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "ifEntry.H"
#include "Switch.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(ifEntry, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        ifEntry,
        execute,
        dictionaryIstream,
        if
    );
} // End namespace functionEntries
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::functionEntries::ifEntry::isTrue(ITstream& its)
{
    Switch logic;

    if (its.front().isScalar())
    {
        // Use default rounding tolerance
        logic = Switch(its.front().scalarToken());
    }
    else
    {
        its >> logic;
    }

    return logic;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionEntries::ifEntry::execute
(
    DynamicList<filePos>& stack,
    dictionary& parentDict,
    Istream& is
)
{
    const label nNested = stack.size();

    stack.push_back(filePos(is.name(), is.lineNumber()));

    // Read line
    string line;
    dynamic_cast<ISstream&>(is).getLine(line);
    line += ';';

    IStringStream lineStream(line);
    const primitiveEntry e("ifEntry", parentDict, lineStream);

    const bool doIf = ifEntry::isTrue(e.stream());

    // Info<< "Using #" << typeName << " " << Switch::name(doIf)
    //     << " at line " << stack.back().second()
    //     << " in file " << stack.back().first() << endl;

    const bool ok = ifeqEntry::execute(doIf, stack, parentDict, is);

    if (stack.size() != nNested)
    {
        FatalIOErrorInFunction(parentDict)
            << "Did not find matching #endif for condition starting"
            << " at line " << stack.back().second()
            << " in file " << stack.back().first() << exit(FatalIOError);
    }

    return ok;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::ifEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    DynamicList<filePos> stack(10);
    return execute(stack, parentDict, is);
}


// ************************************************************************* //
