/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

#include "faMesh.H"
#include "polyMesh.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

bool Foam::faMesh::hasSystemFiles
(
    const word& meshName,
    const polyMesh& pMesh
)
{
    // Expect
    // - system/finite-area/<region>/faSchemes
    // - system/finite-area/<region>/faSolution

    // The directory relative to polyMesh (not Time)
    const fileName relativeDir
    (
        faMesh::prefix() / polyMesh::regionName(meshName)
    );

    DebugInfo<< "check system files: " << relativeDir << nl;

    IOobject systemIOobject
    (
        "any-name",
        pMesh.time().system(),
        relativeDir,
        pMesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    const fileOperation& fp = Foam::fileHandler();

    bool looksValid = true;

    // Global files: system/{faSchemes,faSolution}
    for
    (
        const word& expect
      : List<word>
        ({
            {"faSchemes"},
            {"faSolution"}
        })
    )
    {
        systemIOobject.resetHeader(expect);

        fileName found
        (
            fp.filePath
            (
                true,  // global
                systemIOobject,
                expect // typeName (ununsed?)
            )
        );

        if (found.empty())
        {
            looksValid = false;
        }
    }

    // Only needed on master
    Pstream::broadcast(looksValid);

    return looksValid;
}


bool Foam::faMesh::hasMeshFiles
(
    const word& meshName,
    const polyMesh& pMesh
)
{
    // As well as system/finite-area/{faSchemes,faSolution}
    //
    // expect these:
    // - instance/finite-area/<region>/faMesh/faceLabels
    // - instance/finite-area/<region>/faMesh/faBoundary


    // Not required...
    // bool looksValid = hasSystemFiles(meshName, pMesh);

    bool looksValid = true;

    // The mesh directory relative to polyMesh (not Time)
    const fileName relativeDir
    (
        faMesh::meshDir(word::null, meshName)
    );

    if (looksValid)
    {
        DebugInfo<< "check mesh files: " << relativeDir << nl;

        const fileOperation& fp = Foam::fileHandler();

        // The geometry instance for faMesh/faceLabels
        // Must use READ_IF_PRESENT to avoid aborting if not available

        const word instance = pMesh.time().findInstance
        (
            // Searching from Time, so need polyMesh region too
            pMesh.regionName()/relativeDir,
            "faceLabels",
            IOobject::READ_IF_PRESENT
        );

        IOobject meshIOobject
        (
            "any-name",
            instance,
            relativeDir,
            pMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        for
        (
            const wordPair& expect
          : List<wordPair>
            ({
                {"faceLabels", "labelList"},
                {"faBoundary", "faBoundaryMesh"}
            })
        )
        {
            const word& dataFile = expect.first();
            const word& dataClass = expect.second();

            meshIOobject.resetHeader(dataFile);

            fileName found
            (
                fp.filePath
                (
                    false,      // non-global
                    meshIOobject,
                    dataClass   // typeName (ununsed?)
                )
            );

            if (found.empty())
            {
                looksValid = false;
            }
        }

        // Everybody needs it, or they all fail
        Pstream::reduceAnd(looksValid);
    }

    return looksValid;
}


Foam::autoPtr<Foam::faMesh> Foam::faMesh::TryNew
(
    const word& meshName,
    const polyMesh& pMesh
)
{
    if (faMesh::hasMeshFiles(meshName, pMesh))
    {
        return autoPtr<faMesh>::New(meshName, pMesh);
    }

    return nullptr;
}


Foam::autoPtr<Foam::faMesh> Foam::faMesh::TryNew
(
    const polyMesh& pMesh
)
{
    return TryNew(polyMesh::defaultRegion, pMesh);
}


// ************************************************************************* //
