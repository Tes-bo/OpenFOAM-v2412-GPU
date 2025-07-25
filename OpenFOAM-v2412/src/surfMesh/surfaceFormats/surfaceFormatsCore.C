/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "surfaceFormatsCore.H"
#include "Time.H"
#include "ListOps.H"
#include "surfMesh.H"
#include "stringListOps.H"   // For stringListOps::findMatching()

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::fileFormats::surfaceFormatsCore::nativeExt("ofs");


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::string Foam::fileFormats::surfaceFormatsCore::getLineNoComment
(
    ISstream& is,
    const char comment
)
{
    Foam::string line;
    do
    {
        is.getLine(line);
    }
    while ((line.empty() || line[0] == comment) && is.good());

    return line;
}


Foam::labelList Foam::fileFormats::surfaceFormatsCore::getSelectedPatches
(
    const surfZoneList& patches,
    const wordRes& allow,
    const wordRes& deny
)
{
    return
        stringListOps::findMatching
        (
            patches,
            allow,
            deny,
            nameOp<surfZone>()
        );
}


#if 0
Foam::fileName Foam::fileFormats::surfaceFormatsCore::localMeshFileName
(
    const word& surfName
)
{
    const word name(surfName.size() ? surfName : surfaceRegistry::defaultName);

    return fileName
    (
        surfaceRegistry::prefix/name/surfMesh::meshSubDir
      / name + "." + nativeExt
    );
}


Foam::fileName Foam::fileFormats::surfaceFormatsCore::findMeshInstance
(
    const Time& t,
    const word& surfName
)
{
    fileName localName = localMeshFileName(surfName);

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = t.times();
    label instanceI;

    for (instanceI = ts.size()-1; instanceI >= 0; --instanceI)
    {
        if (ts[instanceI].value() <= t.timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (instanceI >= 0)
    {
        for (label i = instanceI; i >= 0; --i)
        {
            if (isFile(t.path()/ts[i].name()/localName))
            {
                return ts[i].name();
            }
        }
    }

    return t.constant();
}


Foam::fileName Foam::fileFormats::surfaceFormatsCore::findMeshFile
(
    const Time& t,
    const word& surfName
)
{
    fileName localName = localMeshFileName(surfName);

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = t.times();
    label instanceI;

    for (instanceI = ts.size()-1; instanceI >= 0; --instanceI)
    {
        if (ts[instanceI].value() <= t.timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (instanceI >= 0)
    {
        for (label i = instanceI; i >= 0; --i)
        {
            fileName testName(t.path()/ts[i].name()/localName);

            if (isFile(testName))
            {
                return testName;
            }
        }
    }

    // fallback to "constant"
    return t.path()/t.constant()/localName;
}
#endif


Foam::fileName Foam::fileFormats::surfaceFormatsCore::relativeFilePath
(
    const IOobject& io,
    const fileName& f,
    const bool isGlobal
)
{
    fileName fName(f);
    fName.expand();
    if (!fName.isAbsolute())
    {
        // Is the specified file:
        // - local to the cwd?
        // - local to the case dir?
        // - or just another name?
        fName = fileHandler().filePath
        (
            isGlobal,
            IOobject(io, fName),
            word::null
        );
    }
    return fName;
}


Foam::fileName Foam::fileFormats::surfaceFormatsCore::findFile
(
    const IOobject& io,
    const bool isGlobal
)
{
    fileName fName
    (
        isGlobal
      ? io.globalFilePath(word::null)
      : io.localFilePath(word::null)
    );

    if (!exists(fName))
    {
        fName.clear();
    }

    return fName;
}


Foam::fileName Foam::fileFormats::surfaceFormatsCore::findFile
(
    const IOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
{
    fileName fName;
    if (dict.readIfPresent("file", fName, keyType::LITERAL))
    {
        fName = relativeFilePath(io, fName, isGlobal);
    }
    else
    {
        fName =
        (
            isGlobal
          ? io.globalFilePath(word::null)
          : io.localFilePath(word::null)
        );
    }

    if (!exists(fName))
    {
        fName.clear();
    }

    return fName;
}


Foam::fileName Foam::fileFormats::surfaceFormatsCore::checkFile
(
    const IOobject& io,
    const bool isGlobal
)
{
    fileName fName
    (
        isGlobal
      ? io.globalFilePath(word::null)
      : io.localFilePath(word::null)
    );

    if (fName.empty())
    {
        FatalErrorInFunction
            << "Cannot find surface starting from "
            << io.objectPath() << nl
            << exit(FatalError);
    }

    return fName;
}


Foam::fileName Foam::fileFormats::surfaceFormatsCore::checkFile
(
    const IOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
{
    fileName fName;
    if (dict.readIfPresent("file", fName, keyType::LITERAL))
    {
        const fileName rawFName(fName);

        fName = relativeFilePath(io, rawFName, isGlobal);

        if (!exists(fName))
        {
            FatalErrorInFunction
                << "Cannot find surface " << rawFName
                << " starting from " << io.objectPath() << nl
                << exit(FatalError);
        }
    }
    else
    {
        fName =
        (
            isGlobal
          ? io.globalFilePath(word::null)
          : io.localFilePath(word::null)
        );

        if (!exists(fName))
        {
            FatalErrorInFunction
                << "Cannot find surface starting from "
                << io.objectPath() << nl
                << exit(FatalError);
        }
    }

    return fName;
}


bool Foam::fileFormats::surfaceFormatsCore::checkSupport
(
    const wordHashSet& available,
    const word& fileType,
    const bool verbose,
    const char* functionName
)
{
    if (available.found(fileType))
    {
        return true;
    }
    else if (verbose)
    {
        Info<< "Unknown file type";

        if (functionName)
        {
            Info<< " for " << functionName;
        }

        Info<< " : " << fileType << nl
            << "Valid types: " << flatOutput(available.sortedToc()) << nl
            << nl;
    }

    return false;
}


// ************************************************************************* //
