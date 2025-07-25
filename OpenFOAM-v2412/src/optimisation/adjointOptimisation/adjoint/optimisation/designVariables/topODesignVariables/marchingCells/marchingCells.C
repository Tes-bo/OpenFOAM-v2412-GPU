/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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

#include "marchingCells.H"
#include "fvMeshSubset.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(marchingCells, 1);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::marchingCells::initialise()
{
    label nChangedFaces(0);
    labelList changedFaces(mesh_.nFaces(), -1);

    // Gather initial faces from the seeding patches
    for (const label patchI : seedPatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const label start = patch.start();
        forAll(patch, faceI)
        {
            changedFaces[nChangedFaces++] = start + faceI;
        }
    }
    // Gather the boundary faces from the seeding cell zones
    for (label cellZoneID : seedCellZoneIDs_)
    {
        const labelList& zoneCells = mesh_.cellZones()[cellZoneID];
        // Create sub mesh from the given zone
        fvMeshSubset subSetMesh(mesh_, zoneCells);
        const fvMesh& subMesh = subSetMesh.subMesh();
        // Get the exposed faces of the subMesh and use them as seeds to
        // patchWave
        const word patchName(fvMeshSubset::exposedPatchName);
        const label patchID(subMesh.boundaryMesh().findPatchID(patchName));
        const polyPatch& subMeshPatch = subMesh.boundaryMesh()[patchID];
        // Map to faces of the original mesh
        const labelList& faceMap = subSetMesh.faceMap();
        const label start = subMeshPatch.start();
        // Go through a primitivePatch field (faceCentres) since the type
        // of the patch containing the exposed faces is empty and a zero size
        // is returned for the fields accessed through fvPatch
        forAll(subMeshPatch.faceCentres(), faceI)
        {
            changedFaces[nChangedFaces++] = faceMap[start + faceI];
        }
    }
    // Gather faces from the seeding face zones
    for (label faceZoneID : seedFaceZoneIDs_)
    {
        const labelList& zoneFaces = mesh_.faceZones()[faceZoneID];
        for (label faceI : zoneFaces)
        {
            changedFaces[nChangedFaces++] = faceI;
        }
    }
    changedFaces.setSize(nChangedFaces);
    List<wallPointData<bool>> changedFacesInfo(nChangedFaces);
    const vectorField& faceCentres = mesh_.faceCentres();
    forAll(changedFaces, fI)
    {
        changedFacesInfo[fI] =
            wallPointData<bool>(faceCentres[changedFaces[fI]], true, 0.0);
    }

    meshWave_.setFaceInfo(changedFaces, changedFacesInfo);

    // Initialization has been completed
    initialised_ = true;
}


void Foam::marchingCells::appendSeedCell(const label cellID)
{
    if (!isFixedCell_[cellID])
    {
        isActiveCell_[cellID] = true;
        activeCells_.append(cellID);
    }
}


void Foam::marchingCells::march
(
    label nVisited,
    const label cI,
    labelList& newlyAddedCells
)
{
    ++nVisited;
    if (nVisited < marchingStep_ + 1)
    {
        const labelList& cellCells = mesh_.cellCells()[cI];
        for (label neiCellI : cellCells)
        {
            if (!isFixedCell_[neiCellI])
            {
                if (!isActiveCell_[neiCellI])
                {
                    isActiveCell_[neiCellI] = true;
                    newlyAddedCells.append(neiCellI);
                }
                march(nVisited, neiCellI, newlyAddedCells);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::marchingCells::marchingCells
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    seedPatches_
    (
        mesh_.boundaryMesh().patchSet
        (
            dict.getOrDefault<wordRes>("seedPatches", wordRes(0))
        )
    ),
    seedCellZoneIDs_
    (
        mesh_.cellZones().indices
        (
            dict.getOrDefault<wordRes>("seedCellZones", wordRes(0))
        )
    ),
    seedFaceZoneIDs_
    (
        mesh_.faceZones().indices
        (
            dict.getOrDefault<wordRes>("seedFaceZones", wordRes(0))
        )
    ),
    marchingStep_(dict.get<label>("marchingStep")),
    isActiveCell_(mesh_.nCells(), false),
    isFixedCell_(mesh_.nCells(), false),
    activeCells_(0),
    addedCells_(0),
    initialised_(false),
    nIters_(0),
    allFaceInfo_(mesh_.nFaces()),
    allCellInfo_(mesh.nCells()),
    meshWave_(mesh_, allFaceInfo_, allCellInfo_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::marchingCells::update(const label iters)
{
    // If the first seeding has been performed, do it now
    if (!initialised_)
    {
        initialise();
    }

    // Iterate the meshWave algorithm
    meshWave_.iterate(iters*marchingStep_);

    // Grab the newly added cells
    addedCells_ = labelList(mesh_.nCells(), -1);
    label nAddedCells(0);
    forAll(allCellInfo_, cI)
    {
        if (allCellInfo_[cI].data() && !isActiveCell_[cI] && !isFixedCell_[cI])
        {
            addedCells_[nAddedCells++] = cI;
            isActiveCell_[cI] = true;
        }
    }
    addedCells_.setSize(nAddedCells);

    // Add newly found cells to the list of activeCells
    activeCells_.append(addedCells_);

    // Write cell set
    if (debug)
    {
        cellSet activeCellList
        (
            mesh_,
            "activeCells" + name(nIters_),
            activeCells_.size()
        );
        for (label cellI : activeCells_)
        {
            activeCellList.insert(cellI);
        }
        activeCellList.write();
    }

    nIters_ += iters;
}


void Foam::marchingCells::addFixedCells
(
    const cellZoneMesh& cellZoneMesh,
    const labelList& fixedCellZoneIDs
)
{
    for (const label cellZoneID : fixedCellZoneIDs)
    {
        for (const label cI : cellZoneMesh[cellZoneID])
        {
            isFixedCell_[cI] = true;
            isActiveCell_[cI] = false;
        }
    }
}


void Foam::marchingCells::addFixedCells(const labelList& fixedCellIDs)
{
    for (const label cI : fixedCellIDs)
    {
        isFixedCell_[cI] = true;
        isActiveCell_[cI] = false;
    }
}


// ************************************************************************* //
