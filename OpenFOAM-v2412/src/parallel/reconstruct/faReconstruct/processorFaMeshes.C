/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "processorFaMeshes.H"
#include "Time.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorFaMeshes::read()
{
    // Make sure to clear (and hence unregister) any previously loaded meshes
    // and fields

    boundaryProcAddressing_.free();
    faceProcAddressing_.free();
    edgeProcAddressing_.free();
    pointProcAddressing_.free();
    meshes_.free();


    forAll(fvMeshes_, proci)
    {
        // Construct from polyMesh IO information
        meshes_.emplace_set(proci, fvMeshes_[proci]);

        // Read the addressing information

        IOobject ioAddr
        (
            "procAddressing",
            "constant",  // Placeholder
            faMesh::meshSubDir,
            meshes_[proci].thisDb(),
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE
        );

        const auto& runTime = meshes_[proci].thisDb().time();
        const auto& meshDir = meshes_[proci].meshDir();

        // pointProcAddressing (faMesh)
        ioAddr.rename("pointProcAddressing");
        ioAddr.instance() = runTime.findInstance(meshDir, ioAddr.name());
        pointProcAddressing_.emplace_set(proci, ioAddr);

        // edgeProcAddressing (faMesh)
        ioAddr.rename("edgeProcAddressing");
        ioAddr.instance() = runTime.findInstance(meshDir, ioAddr.name());
        edgeProcAddressing_.emplace_set(proci, ioAddr);

        // faceProcAddressing (faMesh)
        ioAddr.rename("faceProcAddressing");
        ioAddr.instance() = runTime.findInstance(meshDir, ioAddr.name());
        faceProcAddressing_.emplace_set(proci, ioAddr);

        // boundaryProcAddressing (faMesh)
        ioAddr.rename("boundaryProcAddressing");
        ioAddr.instance() = runTime.findInstance(meshDir, ioAddr.name());
        boundaryProcAddressing_.emplace_set(proci, ioAddr);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorFaMeshes::processorFaMeshes
(
    const UPtrList<fvMesh>& procFvMeshes
)
:
    fvMeshes_(procFvMeshes),
    meshes_(procFvMeshes.size()),
    pointProcAddressing_(meshes_.size()),
    edgeProcAddressing_(meshes_.size()),
    faceProcAddressing_(meshes_.size()),
    boundaryProcAddressing_(meshes_.size())
{
    read();
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::processorFaMeshes::removeFiles(const faMesh& mesh)
{
    IOobject io
    (
        "procAddressing",
        mesh.facesInstance(),
        faMesh::meshSubDir,
        mesh.thisDb()
    );

    // procAddressing
    fileHandler().rm(fileHandler().filePath(io.objectPath()));

    // pointProcAddressing
    io.rename("pointProcAddressing");
    fileHandler().rm(fileHandler().filePath(io.objectPath()));

    // edgeProcAddressing
    io.rename("edgeProcAddressing");
    fileHandler().rm(fileHandler().filePath(io.objectPath()));

    // faceProcAddressing
    io.rename("faceProcAddressing");
    fileHandler().rm(fileHandler().filePath(io.objectPath()));

    // boundaryProcAddressing
    io.rename("boundaryProcAddressing");
    fileHandler().rm(fileHandler().filePath(io.objectPath()));
}


// ************************************************************************* //
