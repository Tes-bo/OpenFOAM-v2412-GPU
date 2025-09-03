/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "processorMeshes.H"
#include "Time.H"
#include "IndirectList.H"
#include "primitiveMesh.H"
#include "OSspecific.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorMeshes, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorMeshes::read()
{
    // Make sure to clear (and hence unregister) any previously loaded meshes
    // and fields
    pBoundaryProcAddressing_.free();
    pMeshes_.free();

    boundaryProcAddressing_.free();
    cellProcAddressing_.free();
    faceProcAddressing_.free();
    pointProcAddressing_.free();
    meshes_.free();

    forAll(databases_, proci)
    {
        meshes_.emplace_set
        (
            proci,
            IOobject
            (
                meshName_,
                databases_[proci].timeName(),
                databases_[proci]
            )
        );

        // Read the addressing information

        IOobject ioAddr
        (
            "procAddressing",
            meshes_[proci].facesInstance(),
            polyMesh::meshSubDir,
            meshes_[proci].thisDb(),
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE
        );

        // pointProcAddressing (polyMesh)
        ioAddr.rename("pointProcAddressing");
        pointProcAddressing_.emplace_set(proci, ioAddr);

        // faceProcAddressing (polyMesh)
        ioAddr.rename("faceProcAddressing");
        faceProcAddressing_.emplace_set(proci, ioAddr);

        // cellProcAddressing (polyMesh)
        ioAddr.rename("cellProcAddressing");
        cellProcAddressing_.emplace_set(proci, ioAddr);

        // boundaryProcAddressing (polyMesh)
        ioAddr.rename("boundaryProcAddressing");
        boundaryProcAddressing_.emplace_set(proci, ioAddr);


        // pointMesh
        // ~~~~~~~~~

        pMeshes_.set
        (
            proci,
            new pointMesh
            (
                meshes_[proci],
                IOobject::READ_IF_PRESENT
            )
        );

        pBoundaryProcAddressing_.set
        (
            proci,
            autoPtr<labelIOList>::New
            (
                IOobject
                (
                    "boundaryProcAddressing",
                    meshes_[proci].facesInstance(),
                    polyMesh::meshSubDir/pointMesh::meshSubDir,
                    pMeshes_[proci].thisDb(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                boundaryProcAddressing_[proci]
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorMeshes::processorMeshes
(
    PtrList<Time>& databases,
    const word& meshName
)
:
    meshName_(meshName),
    databases_(databases),
    meshes_(databases.size()),
    pointProcAddressing_(databases.size()),
    faceProcAddressing_(databases.size()),
    cellProcAddressing_(databases.size()),
    boundaryProcAddressing_(databases.size()),
    pMeshes_(databases.size()),
    pBoundaryProcAddressing_(databases.size())
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::polyMesh::readUpdateState Foam::processorMeshes::readUpdate()
{
    polyMesh::readUpdateState stat = polyMesh::UNCHANGED;

    forAll(databases_, proci)
    {
        // Check if any new meshes need to be read.
        polyMesh::readUpdateState procStat = meshes_[proci].readUpdate();

        /*
        if (procStat != polyMesh::UNCHANGED)
        {
            Info<< "Processor " << proci
                << " at time " << databases_[proci].timeName()
                << " detected mesh change " << procStat
                << endl;
        }
        */

        // Combine into overall mesh change status
        if (stat == polyMesh::UNCHANGED)
        {
            stat = procStat;
        }
        else if (stat != procStat)
        {
            FatalErrorInFunction
                << "Processor " << proci
                << " has a different polyMesh at time "
                << databases_[proci].timeName()
                << " compared to any previous processors." << nl
                << "Please check time " << databases_[proci].timeName()
                << " directories on all processors for consistent"
                << " mesh files."
                << exit(FatalError);
        }
    }

    if
    (
        stat == polyMesh::TOPO_CHANGE
     || stat == polyMesh::TOPO_PATCH_CHANGE
    )
    {
        // Reread all meshes and addressing
        read();
    }
    return stat;
}


void Foam::processorMeshes::reconstructPoints(fvMesh& mesh)
{
    // Read the field for all the processors
    PtrList<pointIOField> procsPoints(meshes_.size());

    forAll(meshes_, proci)
    {
        procsPoints.set
        (
            proci,
            new pointIOField
            (
                IOobject
                (
                    "points",
                    meshes_[proci].time().timeName(),
                    polyMesh::meshSubDir,
                    meshes_[proci].thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            )
        );
    }

    // Create the new points
    vectorField newPoints(mesh.nPoints());

    forAll(meshes_, proci)
    {
        const vectorField& procPoints = procsPoints[proci];

        const labelList& pointProcAddr = pointProcAddressing_[proci];

        if (pointProcAddr.size() != procPoints.size())
        {
            FatalErrorInFunction
                << "problem :"
                << " pointProcAddr:" << pointProcAddr.size()
                << " procPoints:" << procPoints.size()
                << abort(FatalError);
        }

        UIndirectList<point>(newPoints, pointProcAddr) = procPoints;
        // or: newPoints.rmap(procPoints, pointProcAddr)
    }

    mesh.movePoints(newPoints);
    mesh.write();
}


void Foam::processorMeshes::removeFiles(const polyMesh& mesh)
{
    IOobject io
    (
        "procAddressing",
        mesh.facesInstance(),
        polyMesh::meshSubDir,
        mesh.thisDb()
    );

    // procAddressing
    fileHandler().rm(fileHandler().filePath(io.objectPath()));

    // pointProcAddressing
    io.rename("pointProcAddressing");
    fileHandler().rm(fileHandler().filePath(io.objectPath()));

    // faceProcAddressing
    io.rename("faceProcAddressing");
    fileHandler().rm(fileHandler().filePath(io.objectPath()));

    // cellProcAddressing
    io.rename("cellProcAddressing");
    fileHandler().rm(fileHandler().filePath(io.objectPath()));

    // boundaryProcAddressing
    io.rename("boundaryProcAddressing");
    fileHandler().rm(fileHandler().filePath(io.objectPath()));



    // pointMesh
    // ~~~~~~~~~

    IOobject pointIO
    (
        "boundary",
        mesh.facesInstance(),
        polyMesh::meshSubDir/pointMesh::meshSubDir,
        mesh.thisDb()
    );

    // pointMesh/boundary
    fileHandler().rm(fileHandler().filePath(pointIO.objectPath()));

    // boundaryProcAddressing
    io.rename("boundaryProcAddressing");
    fileHandler().rm(fileHandler().filePath(pointIO.objectPath()));
}


// ************************************************************************* //
