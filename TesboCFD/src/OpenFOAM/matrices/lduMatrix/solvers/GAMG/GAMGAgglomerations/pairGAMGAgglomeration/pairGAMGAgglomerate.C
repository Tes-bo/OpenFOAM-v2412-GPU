/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "pairGAMGAgglomeration.H"
#include "lduAddressing.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::pairGAMGAgglomeration::agglomerate
(
    const label nCellsInCoarsestLevel,
    const label startLevel,
    const scalarField& startFaceWeights,
    const bool doProcessorAgglomerate
)
{
    if (nCells_.size() < maxLevels_)
    {
        // See compactLevels. Make space if not enough
        nCells_.resize(maxLevels_);
        restrictAddressing_.resize(maxLevels_);
        nFaces_.resize(maxLevels_);
        faceRestrictAddressing_.resize(maxLevels_);
        faceFlipMap_.resize(maxLevels_);
        nPatchFaces_.resize(maxLevels_);
        patchFaceRestrictAddressing_.resize(maxLevels_);
        meshLevels_.resize(maxLevels_);
        // Have procCommunicator_ always, even if not procAgglomerating.
        // Use value -1 to indicate nothing is proc-agglomerated
        procCommunicator_.resize(maxLevels_ + 1, -1);
        if (processorAgglomerate())
        {
            procAgglomMap_.resize(maxLevels_);
            agglomProcIDs_.resize(maxLevels_);
            procCommunicator_.resize(maxLevels_);
            procCellOffsets_.resize(maxLevels_);
            procFaceMap_.resize(maxLevels_);
            procBoundaryMap_.resize(maxLevels_);
            procBoundaryFaceMap_.resize(maxLevels_);
        }
    }


    // Start geometric agglomeration from the given faceWeights
    scalarField faceWeights = startFaceWeights;

    // Agglomerate until the required number of cells in the coarsest level
    // is reached

    label nPairLevels = 0;
    label nCreatedLevels = startLevel;

    while (nCreatedLevels < maxLevels_ - 1)
    {
        if (!hasMeshLevel(nCreatedLevels))
        {
            FatalErrorInFunction<< "No mesh at nCreatedLevels:"
                << nCreatedLevels
                << exit(FatalError);
        }

        const auto& fineMesh = meshLevel(nCreatedLevels);


        label nCoarseCells = -1;

        tmp<labelField> finalAgglomPtr = agglomerate
        (
            nCoarseCells,
            fineMesh.lduAddr(),
            faceWeights
        );

        if
        (
            continueAgglomerating
            (
                nCellsInCoarsestLevel,
                finalAgglomPtr().size(),
                nCoarseCells,
                fineMesh.comm()
            )
        )
        {
            nCells_[nCreatedLevels] = nCoarseCells;
            restrictAddressing_.set(nCreatedLevels, finalAgglomPtr);
        }
        else
        {
            break;
        }

        // Create coarse mesh
        agglomerateLduAddressing(nCreatedLevels);

        // Agglomerate the faceWeights field for the next level
        {
            scalarField aggFaceWeights
            (
                meshLevels_[nCreatedLevels].upperAddr().size(),
                0.0
            );

            restrictFaceField
            (
                aggFaceWeights,
                faceWeights,
                nCreatedLevels
            );

            faceWeights = std::move(aggFaceWeights);
        }

        if (nPairLevels % mergeLevels_)
        {
            combineLevels(nCreatedLevels);
        }
        else
        {
            nCreatedLevels++;
        }

        nPairLevels++;
    }

    // Shrink the storage of the levels to those created
    compactLevels(nCreatedLevels, doProcessorAgglomerate);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::pairGAMGAgglomeration::agglomerate
(
    label& nCoarseCells,
    const lduAddressing& fineMatrixAddressing,
    const scalarField& faceWeights
)
{
    const label nFineCells = fineMatrixAddressing.size();

    const labelUList& upperAddr = fineMatrixAddressing.upperAddr();
    const labelUList& lowerAddr = fineMatrixAddressing.lowerAddr();

    // For each cell calculate faces
    labelList cellFaces(upperAddr.size() + lowerAddr.size());
    labelList cellFaceOffsets(nFineCells + 1);

    // memory management
    {
        labelList nNbrs(nFineCells, Zero);

        forAll(upperAddr, facei)
        {
            nNbrs[upperAddr[facei]]++;
        }

        forAll(lowerAddr, facei)
        {
            nNbrs[lowerAddr[facei]]++;
        }

        cellFaceOffsets[0] = 0;
        forAll(nNbrs, celli)
        {
            cellFaceOffsets[celli+1] = cellFaceOffsets[celli] + nNbrs[celli];
        }

        // reset the whole list to use as counter
        nNbrs = 0;

        forAll(upperAddr, facei)
        {
            cellFaces
            [
                cellFaceOffsets[upperAddr[facei]] + nNbrs[upperAddr[facei]]
            ] = facei;

            nNbrs[upperAddr[facei]]++;
        }

        forAll(lowerAddr, facei)
        {
            cellFaces
            [
                cellFaceOffsets[lowerAddr[facei]] + nNbrs[lowerAddr[facei]]
            ] = facei;

            nNbrs[lowerAddr[facei]]++;
        }
    }


    // go through the faces and create clusters

    auto tcoarseCellMap = tmp<labelField>::New(nFineCells, -1);
    auto& coarseCellMap = tcoarseCellMap.ref();

    nCoarseCells = 0;
    label celli;

    for (label cellfi=0; cellfi<nFineCells; cellfi++)
    {
        // Change cell ordering depending on direction for this level
        celli = forward_ ? cellfi : nFineCells - cellfi - 1;

        if (coarseCellMap[celli] < 0)
        {
            label matchFaceNo = -1;
            scalar maxFaceWeight = -GREAT;

            // check faces to find ungrouped neighbour with largest face weight
            for
            (
                label faceOs=cellFaceOffsets[celli];
                faceOs<cellFaceOffsets[celli+1];
                faceOs++
            )
            {
                label facei = cellFaces[faceOs];

                // I don't know whether the current cell is owner or neighbour.
                // Therefore I'll check both sides
                if
                (
                    coarseCellMap[upperAddr[facei]] < 0
                 && coarseCellMap[lowerAddr[facei]] < 0
                 && faceWeights[facei] > maxFaceWeight
                )
                {
                    // Match found. Pick up all the necessary data
                    matchFaceNo = facei;
                    maxFaceWeight = faceWeights[facei];
                }
            }

            if (matchFaceNo >= 0)
            {
                // Make a new group
                coarseCellMap[upperAddr[matchFaceNo]] = nCoarseCells;
                coarseCellMap[lowerAddr[matchFaceNo]] = nCoarseCells;
                nCoarseCells++;
            }
            else
            {
                // No match. Find the best neighbouring cluster and
                // put the cell there
                label clusterMatchFaceNo = -1;
                scalar clusterMaxFaceCoeff = -GREAT;

                for
                (
                    label faceOs=cellFaceOffsets[celli];
                    faceOs<cellFaceOffsets[celli+1];
                    faceOs++
                )
                {
                    label facei = cellFaces[faceOs];

                    if (faceWeights[facei] > clusterMaxFaceCoeff)
                    {
                        clusterMatchFaceNo = facei;
                        clusterMaxFaceCoeff = faceWeights[facei];
                    }
                }

                if (clusterMatchFaceNo >= 0)
                {
                    // Add the cell to the best cluster
                    coarseCellMap[celli] = max
                    (
                        coarseCellMap[upperAddr[clusterMatchFaceNo]],
                        coarseCellMap[lowerAddr[clusterMatchFaceNo]]
                    );
                }
            }
        }
    }

    // Check that all cells are part of clusters,
    // if not create single-cell "clusters" for each
    for (label cellfi=0; cellfi<nFineCells; cellfi++)
    {
        // Change cell ordering depending on direction for this level
        celli = forward_ ? cellfi : nFineCells - cellfi - 1;

        if (coarseCellMap[celli] < 0)
        {
            coarseCellMap[celli] = nCoarseCells;
            nCoarseCells++;
        }
    }

    if (!forward_)
    {
        nCoarseCells--;

        forAll(coarseCellMap, celli)
        {
            coarseCellMap[celli] = nCoarseCells - coarseCellMap[celli];
        }

        nCoarseCells++;
    }

    // Reverse the map ordering for the next level
    // to improve the next level of agglomeration
    forward_ = !forward_;

    return tcoarseCellMap;
}


// ************************************************************************* //
