/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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

#include "includeEtcEntry.H"
#include "etcFiles.H"
#include "stringOps.H"
#include "IFstream.H"
#include "IOstreams.H"
#include "UPstream.H"
#include "fileOperation.H"
#include "regIOobject.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEtcEntry::log(false);


namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEtcEntry,
        execute,
        dictionaryIstream,
        includeEtc
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEtcEntry,
        execute,
        primitiveEntryIstream,
        includeEtc
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        sincludeEtcEntry,
        execute,
        dictionaryIstream,
        sincludeEtc
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        sincludeEtcEntry,
        execute,
        primitiveEntryIstream,
        sincludeEtc
    );
} // End namespace functionEntries
} // End namespace Foam


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEtcEntry::resolveEtcFile
(
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

    // Search etc directories for the file
    return Foam::findEtcFile(fName);
}


bool Foam::functionEntries::includeEtcEntry::execute
(
    const bool mandatory,
    dictionary& parentDict,
    Istream& is
)
{
    const regIOobject* rioPtr = isA<regIOobject>(parentDict.topDict());

    const label oldComm
    (
        rioPtr && rioPtr->global()
      ? fileHandler().comm(UPstream::worldComm)
      : fileHandler().comm()
    );

    const fileName rawName(is);
    const fileName fName(resolveEtcFile(rawName, parentDict));

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    auto& ifs = *ifsPtr;

    if (ifs)
    {
        if (Foam::functionEntries::includeEtcEntry::log)
        {
            // Report to stdout which file is included
            Info<< fName << nl;
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
        << "Cannot open etc file "
        << (ifs.name().size() ? ifs.name() : rawName)
        << " while reading dictionary " << parentDict.relativeName()
        << exit(FatalIOError);

    return false;
}


bool Foam::functionEntries::includeEtcEntry::execute
(
    const bool mandatory,
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    const regIOobject* rioPtr = isA<regIOobject>(parentDict.topDict());

    const label oldComm
    (
        rioPtr && rioPtr->global()
      ? fileHandler().comm(UPstream::worldComm)
      : fileHandler().comm()
    );

    const fileName rawName(is);
    const fileName fName(resolveEtcFile(rawName, parentDict));

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    auto& ifs = *ifsPtr;

    if (ifs)
    {
        if (Foam::functionEntries::includeEtcEntry::log)
        {
            // Report to stdout which file is included
            Info<< fName << nl;
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
        << "Cannot open etc file "
        << (ifs.name().size() ? ifs.name() : rawName)
        << " while reading dictionary " << parentDict.relativeName()
        << exit(FatalIOError);

    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEtcEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    return includeEtcEntry::execute(true, parentDict, is);
}


bool Foam::functionEntries::includeEtcEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    return includeEtcEntry::execute(true, parentDict, entry, is);
}


bool Foam::functionEntries::sincludeEtcEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    return includeEtcEntry::execute(false, parentDict, is);
}


bool Foam::functionEntries::sincludeEtcEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    return includeEtcEntry::execute(false, parentDict, entry, is);
}


// ************************************************************************* //
