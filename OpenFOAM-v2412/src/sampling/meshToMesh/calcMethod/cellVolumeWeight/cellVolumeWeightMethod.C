/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015 OpenCFD Ltd.
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

#include "cellVolumeWeightMethod.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellVolumeWeightMethod, 0);
    addToRunTimeSelectionTable
    (
        meshToMeshMethod,
        cellVolumeWeightMethod,
        components
    );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::cellVolumeWeightMethod::findInitialSeeds
(
    const labelList& srcCellIDs,
    const boolList& mapFlag,
    const label startSeedI,
    label& srcSeedI,
    label& tgtSeedI
) const
{
    const cellList& srcCells = src_.cells();
    const faceList& srcFaces = src_.faces();
    const pointField& srcPts = src_.points();

    for (label i = startSeedI; i < srcCellIDs.size(); ++i)
    {
        const label srcI = srcCellIDs[i];

        if (mapFlag[srcI])
        {
            const pointField pts(srcCells[srcI].points(srcFaces, srcPts));

            for (const point& pt : pts)
            {
                const label tgtI = tgt_.cellTree().findInside(pt);

                if (tgtI != -1 && intersect(srcI, tgtI))
                {
                    srcSeedI = srcI;
                    tgtSeedI = tgtI;

                    return true;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "could not find starting seed" << endl;
    }

    return false;
}


void Foam::cellVolumeWeightMethod::calculateAddressing
(
    labelListList& srcToTgtCellAddr,
    scalarListList& srcToTgtCellWght,
    labelListList& tgtToSrcCellAddr,
    scalarListList& tgtToSrcCellWght,
    const label srcSeedI,
    const label tgtSeedI,
    const labelList& srcCellIDs,
    boolList& mapFlag,
    label& startSeedI
)
{
    label srcCelli = srcSeedI;
    label tgtCelli = tgtSeedI;

    List<DynamicList<label>> srcToTgtAddr(src_.nCells());
    List<DynamicList<scalar>> srcToTgtWght(src_.nCells());

    List<DynamicList<label>> tgtToSrcAddr(tgt_.nCells());
    List<DynamicList<scalar>> tgtToSrcWght(tgt_.nCells());

    // List of tgt cell neighbour cells
    DynamicList<label> queuedCells(10);

    // List of tgt cells currently visited for srcCellI to avoid multiple hits
    DynamicList<label> visitedCells(10);

    // list to keep track of tgt cells used to seed src cells
    labelList seedCells(src_.nCells(), -1);
    seedCells[srcCelli] = tgtCelli;

    const scalarField& srcVol = src_.cellVolumes();

    do
    {
        queuedCells.clear();
        visitedCells.clear();

        // Initial target cell and neighbours
        queuedCells.push_back(tgtCelli);
        appendNbrCells(tgtCelli, tgt_, visitedCells, queuedCells);

        while (!queuedCells.empty())
        {
            // Process new target cell as LIFO
            tgtCelli = queuedCells.back();
            queuedCells.pop_back();
            visitedCells.push_back(tgtCelli);

            scalar vol = interVol(srcCelli, tgtCelli);

            // accumulate addressing and weights for valid intersection
            if (vol/srcVol[srcCelli] > tolerance_)
            {
                // store src/tgt cell pair
                srcToTgtAddr[srcCelli].append(tgtCelli);
                srcToTgtWght[srcCelli].append(vol);

                tgtToSrcAddr[tgtCelli].append(srcCelli);
                tgtToSrcWght[tgtCelli].append(vol);

                appendNbrCells(tgtCelli, tgt_, visitedCells, queuedCells);

                // accumulate intersection volume
                V_ += vol;
            }
        }

        mapFlag[srcCelli] = false;

        // find new source seed cell
        setNextCells
        (
            startSeedI,
            srcCelli,
            tgtCelli,
            srcCellIDs,
            mapFlag,
            visitedCells,
            seedCells
        );
    }
    while (srcCelli != -1);

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr, i)
    {
        srcToTgtCellAddr[i].transfer(srcToTgtAddr[i]);
        srcToTgtCellWght[i].transfer(srcToTgtWght[i]);
    }

    forAll(tgtToSrcCellAddr, i)
    {
        tgtToSrcCellAddr[i].transfer(tgtToSrcAddr[i]);
        tgtToSrcCellWght[i].transfer(tgtToSrcWght[i]);
    }


    if (debug%2)
    {
        // At this point the overlaps are still in volume so we could
        // get out the relative error
        forAll(srcToTgtCellAddr, cellI)
        {
            scalar srcVol = src_.cellVolumes()[cellI];
            scalar tgtVol = sum(srcToTgtCellWght[cellI]);

            if (mag(srcVol) > ROOTVSMALL && mag((tgtVol-srcVol)/srcVol) > 1e-6)
            {
                WarningInFunction
                    << "At cell " << cellI << " cc:"
                    << src_.cellCentres()[cellI]
                    << " vol:" << srcVol
                    << " total overlap volume:" << tgtVol
                    << endl;
            }
        }

        forAll(tgtToSrcCellAddr, cellI)
        {
            scalar tgtVol = tgt_.cellVolumes()[cellI];
            scalar srcVol = sum(tgtToSrcCellWght[cellI]);

            if (mag(tgtVol) > ROOTVSMALL && mag((srcVol-tgtVol)/tgtVol) > 1e-6)
            {
                WarningInFunction
                    << "At cell " << cellI << " cc:"
                    << tgt_.cellCentres()[cellI]
                    << " vol:" << tgtVol
                    << " total overlap volume:" << srcVol
                    << endl;
            }
        }
    }
}


void Foam::cellVolumeWeightMethod::setNextCells
(
    label& startSeedI,
    label& srcCelli,
    label& tgtCelli,
    const labelList& srcCellIDs,
    const boolList& mapFlag,
    const DynamicList<label>& visitedCells,
    labelList& seedCells
) const
{
    const labelList& srcNbrCells = src_.cellCells()[srcCelli];

    // set possible seeds for later use by querying all src cell neighbours
    // with all visited target cells
    bool valuesSet = false;
    forAll(srcNbrCells, i)
    {
        label cellS = srcNbrCells[i];

        if (mapFlag[cellS] && seedCells[cellS] == -1)
        {
            forAll(visitedCells, j)
            {
                label cellT = visitedCells[j];

                if (intersect(cellS, cellT))
                {
                    seedCells[cellS] = cellT;

                    if (!valuesSet)
                    {
                        srcCelli = cellS;
                        tgtCelli = cellT;
                        valuesSet = true;
                    }
                }
            }
        }
    }

    // set next src and tgt cells if not set above
    if (valuesSet)
    {
        return;
    }
    else
    {
        // try to use existing seed
        bool foundNextSeed = false;
        for (label i = startSeedI; i < srcCellIDs.size(); i++)
        {
            label cellS = srcCellIDs[i];

            if (mapFlag[cellS])
            {
                if (!foundNextSeed)
                {
                    startSeedI = i;
                    foundNextSeed = true;
                }

                if (seedCells[cellS] != -1)
                {
                    srcCelli = cellS;
                    tgtCelli = seedCells[cellS];

                    return;
                }
            }
        }

        // perform new search to find match
        if (debug)
        {
            Pout<< "Advancing front stalled: searching for new "
                << "target cell" << endl;
        }

        bool restart =
            findInitialSeeds
            (
                srcCellIDs,
                mapFlag,
                startSeedI,
                srcCelli,
                tgtCelli
            );

        if (restart)
        {
            // successfully found new starting seed-pair
            return;
        }
    }

    // if we have got to here, there are no more src/tgt cell intersections
    srcCelli = -1;
    tgtCelli = -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellVolumeWeightMethod::cellVolumeWeightMethod
(
    const polyMesh& src,
    const polyMesh& tgt
)
:
    meshToMeshMethod(src, tgt)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::cellVolumeWeightMethod::~cellVolumeWeightMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellVolumeWeightMethod::calculate
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    pointListList& srcToTgtVec,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght,
    pointListList& tgtToSrcVec
)
{
    bool ok = initialise
    (
        srcToTgtAddr,
        srcToTgtWght,
        tgtToSrcAddr,
        tgtToSrcWght
    );

    if (!ok)
    {
        return;
    }

    // (potentially) participating source mesh cells
    const labelList srcCellIDs(maskCells());

    // list to keep track of whether src cell can be mapped
    boolList mapFlag(src_.nCells(), false);
    boolUIndList(mapFlag, srcCellIDs) = true;

    // find initial point in tgt mesh
    label srcSeedI = -1;
    label tgtSeedI = -1;
    label startSeedI = 0;

    bool startWalk =
        findInitialSeeds
        (
            srcCellIDs,
            mapFlag,
            startSeedI,
            srcSeedI,
            tgtSeedI
        );

    if (startWalk)
    {
        calculateAddressing
        (
            srcToTgtAddr,
            srcToTgtWght,
            tgtToSrcAddr,
            tgtToSrcWght,
            srcSeedI,
            tgtSeedI,
            srcCellIDs,
            mapFlag,
            startSeedI
        );
    }
    else
    {
        // if meshes are collocated, after inflating the source mesh bounding
        // box tgt mesh cells may be transferred, but may still not overlap
        // with the source mesh
        return;
    }
}


// ************************************************************************* //
