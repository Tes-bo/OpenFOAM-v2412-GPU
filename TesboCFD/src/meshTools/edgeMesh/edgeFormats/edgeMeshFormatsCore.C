/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "edgeMeshFormatsCore.H"
#include "Time.H"
#include "ListOps.H"
#include "edgeMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::fileFormats::edgeMeshFormatsCore::nativeExt("eMesh");


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::string Foam::fileFormats::edgeMeshFormatsCore::getLineNoComment
(
    ISstream& is,
    const char comment
)
{
    string line;
    do
    {
        is.getLine(line);
    }
    while ((line.empty() || line[0] == comment) && is.good());

    return line;
}


#if 0
Foam::fileName Foam::fileFormats::edgeMeshFormatsCore::localMeshFileName
(
    const word& meshName
)
{
    const word name(meshName.size() ? meshName : surfaceRegistry::defaultName);

    return fileName
    (
        surfaceRegistry::prefix/name/surfMesh::meshSubDir
      / name + "." + nativeExt
    );
}


Foam::fileName Foam::fileFormats::edgeMeshFormatsCore::findMeshInstance
(
    const Time& t,
    const word& meshName
)
{
    fileName localName = localMeshFileName(meshName);

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


Foam::fileName Foam::fileFormats::edgeMeshFormatsCore::findMeshFile
(
    const Time& t,
    const word& meshName
)
{
    fileName localName = localMeshFileName(meshName);

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


bool Foam::fileFormats::edgeMeshFormatsCore::checkSupport
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
