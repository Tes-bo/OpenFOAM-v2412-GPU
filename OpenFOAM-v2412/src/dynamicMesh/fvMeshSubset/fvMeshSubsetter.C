/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "fvMeshSubsetter.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "removeCells.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Helper: extract cells-to-remove from cells-to-keep
static bitSet invertCellSelection
(
    const label nCells,
    const bitSet& selectedCells
)
{
    // Work on a copy
    bitSet cellsToRemove(selectedCells);

    // Ensure we have the full range
    cellsToRemove.resize(nCells, false);

    // Invert the selection
    cellsToRemove.flip();

    return cellsToRemove;
}


// Helper: extract cells-to-remove from cells-to-keep
static inline bitSet invertCellSelection
(
    const label nCells,
    const label regioni,
    const labelUList& regions
)
{
    return BitSetOps::create
    (
        nCells,
        regioni,
        regions,
        false  // on=false: invert return cells to remove
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshSubsetter::removeCellsImpl
(
    const bitSet& cellsToRemove,
    const labelList& exposedFaces,
    const labelList& patchIDs,
    const bool syncPar
)
{
    // Clear out all existing maps
    clear();

    // Mesh changing engine.
    polyTopoChange meshMod(baseMesh());

    removeCells cellRemover(baseMesh(), syncPar);

    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        patchIDs,
        meshMod
    );

    // Create mesh, return map from old to new mesh.
    autoPtr<fvMesh> newMeshPtr;
    autoPtr<mapPolyMesh> map = meshMod.makeMesh
    (
        newMeshPtr,
        IOobject
        (
            baseMesh().name(),
            baseMesh().time().timeName(),
            baseMesh().time(),
            IOobject::READ_IF_PRESENT,  // read fv* if present
            IOobject::NO_WRITE
        ),
        baseMesh(),
        syncPar
    );

    reset
    (
        std::move(newMeshPtr),
        labelList(map().pointMap()),
        labelList(map().faceMap()),
        labelList(map().cellMap()),
        identity(baseMesh().boundaryMesh().size())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::fvMeshSubsetter::getExposedFaces
(
    const bitSet& selectedCells,
    const bool syncPar
) const
{
    return
        Foam::removeCells(baseMesh(), syncPar).getExposedFaces
        (
            invertCellSelection(baseMesh().nCells(), selectedCells)
        );
}


Foam::labelList Foam::fvMeshSubsetter::getExposedFaces
(
    const label regioni,
    const labelUList& regions,
    const bool syncPar
) const
{
    return
        Foam::removeCells(baseMesh(), syncPar).getExposedFaces
        (
            invertCellSelection(baseMesh().nCells(), regioni, regions)
        );
}


void Foam::fvMeshSubsetter::setCellSubset
(
    const bitSet& selectedCells,
    const labelList& exposedFaces,
    const labelList& patchIDs,
    const bool syncPar
)
{
    removeCellsImpl
    (
        invertCellSelection(baseMesh().nCells(), selectedCells),
        exposedFaces,
        patchIDs,
        syncPar
    );
}


void Foam::fvMeshSubsetter::setCellSubset
(
    const label regioni,
    const labelList& regions,
    const labelList& exposedFaces,
    const labelList& patchIDs,
    const bool syncCouples
)
{
    removeCellsImpl
    (
        invertCellSelection(baseMesh().nCells(), regioni, regions),
        exposedFaces,
        patchIDs,
        syncCouples
    );
}


// ************************************************************************* //
