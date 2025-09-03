/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "functionEntry.H"
#include "IOstreams.H"
#include "ISstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineMemberFunctionSelectionTable
    (
        functionEntry,
        execute,
        dictionaryIstream
    );

    defineMemberFunctionSelectionTable
    (
        functionEntry,
        execute,
        primitiveEntryIstream
    );

} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::token Foam::functionEntry::readLine(const word& key, Istream& is)
{
    string s;
    dynamic_cast<ISstream&>(is).getLine(s);

    return token(string(key+s), is.lineNumber());
}


bool Foam::functionEntry::continueReadUntilRightBrace
(
    Istream& is,
    std::string& str,
    const bool stripComments
)
{
    return
        dynamic_cast<ISstream&>(is)
            .continueReadUntilRightBrace(str, stripComments);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntry::functionEntry
(
    const word& key,
    const dictionary& dict,
    Istream& is
)
:
    primitiveEntry
    (
        word(key + dict.name() + Foam::name(is.lineNumber())),
        readLine(key, is)
    )
{}


// * * * * * * * * * * * * Member Function Selectors * * * * * * * * * * * * //

bool Foam::functionEntry::execute
(
    const word& functionName,
    dictionary& parentDict,
    Istream& is
)
{
    is.fatalCheck(FUNCTION_NAME);

    if (!executedictionaryIstreamMemberFunctionTablePtr_)
    {
        std::cerr
            << FUNCTION_NAME << nl
            << "Not yet initialized, function = "
            << functionName.c_str() << std::endl;

        // Return true to keep reading
        return true;
    }

    auto* mfuncPtr =
        executedictionaryIstreamMemberFunctionTable(functionName);

    if (!mfuncPtr)
    {
        FatalErrorInFunction
            << "Unknown functionEntry '" << functionName
            << "' in " << is.relativeName()
            << " near line " << is.lineNumber() << nl << nl
            << "Valid functionEntries :" << nl
            << executedictionaryIstreamMemberFunctionTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return mfuncPtr(parentDict, is);
}


bool Foam::functionEntry::execute
(
    const word& functionName,
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    is.fatalCheck(FUNCTION_NAME);

    if (!executeprimitiveEntryIstreamMemberFunctionTablePtr_)
    {
        std::cerr
            << FUNCTION_NAME << nl
            << "Not yet initialized, function = "
            << functionName.c_str() << std::endl;

        // Return true to keep reading anyhow
        return true;
    }

    auto* mfuncPtr =
        executeprimitiveEntryIstreamMemberFunctionTable(functionName);

    if (!mfuncPtr)
    {
        FatalErrorInFunction
            << "Unknown functionEntry '" << functionName
            << "' in " << is.relativeName()
            << " near line " << is.lineNumber() << nl << nl
            << "Valid functionEntries :" << nl
            << executeprimitiveEntryIstreamMemberFunctionTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return mfuncPtr(parentDict, entry, is);
}


void Foam::functionEntry::write(Ostream& os) const
{
    // Contents should be single string token
    const token& tok = operator[](0);
    const string& s = tok.stringToken();

    // Write character-wise for literal output
    for (size_t i = 0; i < s.size(); ++i)
    {
        os.write(s[i]);
    }

    os << nl;
}


// ************************************************************************* //
