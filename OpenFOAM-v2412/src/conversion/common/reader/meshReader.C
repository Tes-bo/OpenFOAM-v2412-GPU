/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "meshReader.H"
#include "Time.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshReader::addCellZones(polyMesh& mesh) const
{
    cellTable_.addCellZones(mesh, cellTableId_);
    warnDuplicates("cellZones", mesh.cellZones().names());
}


void Foam::meshReader::addFaceZones(polyMesh& mesh) const
{
    label nZone = monitoringSets_.size();
    mesh.faceZones().setSize(nZone);

    if (!nZone)
    {
        return;
    }

    nZone = 0;
    forAllConstIters(monitoringSets_, iter)
    {
        Info<< "faceZone " << nZone
            << " (size: " << iter().size() << ") name: "
            << iter.key() << endl;

        mesh.faceZones().set
        (
            nZone,
            new faceZone
            (
                iter.key(),
                iter(),
                false, // none are flipped
                nZone,
                mesh.faceZones()
            )
        );

        nZone++;
    }
    mesh.faceZones().writeOpt(IOobject::AUTO_WRITE);
    warnDuplicates("faceZones", mesh.faceZones().names());
}


Foam::autoPtr<Foam::polyMesh> Foam::meshReader::mesh
(
    const objectRegistry& registry
)
{
    readGeometry();

    Info<< "Creating a polyMesh" << endl;
    createPolyCells();

    Info<< "Number of internal faces: " << nInternalFaces_ << endl;

    createPolyBoundary();
    clearExtraStorage();

    auto meshPtr = autoPtr<polyMesh>::New
    (
        IOobject
        (
            polyMesh::defaultRegion,
            registry.time().constant(),
            registry,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        std::move(points_),
        std::move(meshFaces_),
        std::move(cellPolys_)
    );
    polyMesh& mesh = *meshPtr;

    polyPatchList newPatches(polyBoundaryPatches(mesh));

    // Adding patches also checks the mesh
    mesh.addPatches(newPatches);

    warnDuplicates("boundaries", mesh.boundaryMesh().names());

    addCellZones(mesh);
    addFaceZones(mesh);

    return meshPtr;
}


void Foam::meshReader::writeMesh
(
    const polyMesh& mesh,
    IOstreamOption streamOpt
) const
{
    mesh.removeFiles();

    Info<< "Writing polyMesh" << endl;
    mesh.writeObject(streamOpt, true);
    writeAux(mesh);
}


void Foam::meshReader::clearExtraStorage()
{
    cellFaces_.clear();
    baffleFaces_.clear();
    boundaryIds_.clear();
    baffleIds_.clear();

    pointCellsPtr_.reset(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshReader::meshReader
(
    const fileName& fileOrPrefix,
    const scalar scaleFactor
)
:
    pointCellsPtr_(nullptr),
    interfaces_(0),
    baffleIds_(0),
    cellPolys_(0),
    monitoringSets_(),
    // protected
    geometryFile_(fileOrPrefix),
    scaleFactor_(scaleFactor),
    points_(0),
    origCellId_(0),
    boundaryIds_(0),
    patchTypes_(0),
    patchNames_(0),
    patchPhysicalTypes_(0),
    patchStarts_(0),
    patchSizes_(0),
    nInternalFaces_(0),
    meshFaces_(0),
    cellFaces_(0),
    baffleFaces_(0),
    cellTableId_(0),
    cellTable_()
{
    // Sanity
    if (scaleFactor_ <= VSMALL)
    {
        scaleFactor_ = 1;
    }
}


// ************************************************************************* //
