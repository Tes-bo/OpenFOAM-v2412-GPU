/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "calcEntry.H"
#include "codeStream.H"
#include "dictionary.H"
#include "dynamicCode.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        calcEntry,
        execute,
        dictionaryIstream,
        calc
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        calcEntry,
        execute,
        primitiveEntryIstream,
        calc
    );
} // End namespace functionEntries
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::string Foam::functionEntries::calcEntry::evaluate
(
    const dictionary& parentDict,
    Istream& is
)
{
    DetailInfo
        << "Using #calc - line "
        << is.lineNumber() << " in file "
        << parentDict.relativeName() << nl;

    dynamicCode::checkSecurity
    (
        "functionEntries::calcEntry::evaluate(..)",
        parentDict
    );

    // Read string
    string s(is);

    // Construct codeDict for codeStream
    dictionary codeSubDict;
    codeSubDict.add("code", "os << (" + s + ");");
    dictionary codeDict(parentDict, codeSubDict);

    // Use function to write stream
    OStringStream buf(is.format());
    buf.precision(16);      // Some reasonably high precision

    streamingFunctionType function = getFunction(parentDict, codeDict);
    (*function)(buf, parentDict);

    // Return evaluated content as string
    return buf.str();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::calcEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    IStringStream result(evaluate(parentDict, is));
    entry.read(parentDict, result);

    return true;
}


bool Foam::functionEntries::calcEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    IStringStream result(evaluate(parentDict, is));
    parentDict.read(result);

    return true;
}


// ************************************************************************* //
