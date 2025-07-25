/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "includeEntry.H"
#include "stringOps.H"
#include "IFstream.H"
#include "IOstreams.H"
#include "UPstream.H"
#include "fileOperation.H"
#include "regIOobject.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEntry::log(false);


namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        dictionaryIstream,
        include
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        primitiveEntryIstream,
        include
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        sincludeEntry,
        execute,
        dictionaryIstream,
        sinclude
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        sincludeEntry,
        execute,
        primitiveEntryIstream,
        sinclude
    );

    // Compat 1712 and earlier
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        sincludeEntry,
        execute,
        dictionaryIstream,
        includeIfPresent
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        sincludeEntry,
        execute,
        primitiveEntryIstream,
        includeIfPresent
    );
} // End namespace functionEntries
} // End namespace Foam


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEntry::resolveFile
(
    const fileName& dir,
    const fileName& f,
    const dictionary& dict
)
{
    fileName fName(f);

    // Substitute dictionary and environment variables.
    // Allow empty substitutions.
    stringOps::inplaceExpand(fName, dict, true, true);

    if (fName.empty() || fName.isAbsolute())
    {
        return fName;
    }

    // Relative name
    return dir/fName;
}


bool Foam::functionEntries::includeEntry::execute
(
    const bool mandatory,
    dictionary& parentDict,
    Istream& is
)
{
    const auto* rioPtr = isA<regIOobject>(parentDict.topDict());

    const label oldComm
    (
        rioPtr && rioPtr->global()
      ? fileHandler().comm(UPstream::worldComm)
      : fileHandler().comm()
    );

    const fileName rawName(is);
    const fileName fName(resolveFile(is.name().path(), rawName, parentDict));

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    auto& ifs = *ifsPtr;

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            // Report to stdout which file is included
            Info<< fName << nl;
        }

        // Add watch on included file
        if (rioPtr)
        {
            const_cast<regIOobject&>(*rioPtr).addWatch(fName);
        }

        parentDict.read(ifs);

        fileHandler().comm(oldComm);
        return true;
    }

    fileHandler().comm(oldComm);

    if (!mandatory)
    {
        return true; // Never fails if optional
    }

    FatalIOErrorInFunction(is)
        << "Cannot open include file "
        << (ifs.name().size() ? ifs.name() : rawName)
        << " while reading dictionary " << parentDict.relativeName()
        << exit(FatalIOError);

    return false;
}


bool Foam::functionEntries::includeEntry::execute
(
    const bool mandatory,
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    const auto* rioPtr = isA<regIOobject>(parentDict.topDict());

    const label oldComm
    (
        rioPtr && rioPtr->global()
      ? fileHandler().comm(UPstream::worldComm)
      : fileHandler().comm()
    );

    const fileName rawName(is);
    const fileName fName(resolveFile(is.name().path(), rawName, parentDict));

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    auto& ifs = *ifsPtr;

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            // Report to stdout which file is included
            Info<< fName << nl;
        }

        // Add watch on included file
        if (rioPtr)
        {
            const_cast<regIOobject&>(*rioPtr).addWatch(fName);
        }

        entry.read(parentDict, ifs);

        fileHandler().comm(oldComm);
        return true;
    }

    fileHandler().comm(oldComm);

    if (!mandatory)
    {
        return true; // Never fails if optional
    }

    FatalIOErrorInFunction(is)
        << "Cannot open include file "
        << (ifs.name().size() ? ifs.name() : rawName)
        << " while reading dictionary " << parentDict.relativeName()
        << exit(FatalIOError);

    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    return includeEntry::execute(true, parentDict, is);
}


bool Foam::functionEntries::includeEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    return includeEntry::execute(true, parentDict, entry, is);
}


bool Foam::functionEntries::sincludeEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    return includeEntry::execute(false, parentDict, is);
}


bool Foam::functionEntries::sincludeEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    return includeEntry::execute(false, parentDict, entry, is);
}


// ************************************************************************* //
