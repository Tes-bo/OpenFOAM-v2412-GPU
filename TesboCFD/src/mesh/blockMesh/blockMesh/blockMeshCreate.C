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

#include "blockMesh.H"
#include "cellModel.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockMesh::createPoints() const
{
    const blockList& blocks = *this;

    const vector scaleTot
    (
        prescaling_.x() * scaling_.x(),
        prescaling_.y() * scaling_.y(),
        prescaling_.z() * scaling_.z()
    );

    if (verbose_)
    {
        Info<< "Creating points with scale " << scaleTot << endl;
    }

    points_.resize(nPoints_);

    forAll(blocks, blocki)
    {
        const pointField& blockPoints = blocks[blocki].points();

        const labelSubList pointAddr
        (
            mergeList_,
            blockPoints.size(),
            blockOffsets_[blocki]
        );

        UIndirectList<point>(points_, pointAddr) = blockPoints;

        if (verbose_)
        {
            Info<< "    Block " << blocki << " cell size :" << nl;

            const label v0 = blocks[blocki].pointLabel(0, 0, 0);

            {
                const label nx = blocks[blocki].density().x();
                const label v1 = blocks[blocki].pointLabel(1, 0, 0);
                const label vn = blocks[blocki].pointLabel(nx, 0, 0);
                const label vn1 = blocks[blocki].pointLabel(nx-1, 0, 0);

                const scalar cwBeg = mag(blockPoints[v1] - blockPoints[v0]);
                const scalar cwEnd = mag(blockPoints[vn] - blockPoints[vn1]);

                Info<< "        i : "
                    << cwBeg*scaleTot.x() << " .. "
                    << cwEnd*scaleTot.x() << nl;
            }

            {
                const label ny = blocks[blocki].density().y();
                const label v1 = blocks[blocki].pointLabel(0, 1, 0);
                const label vn = blocks[blocki].pointLabel(0, ny, 0);
                const label vn1 = blocks[blocki].pointLabel(0, ny-1, 0);

                const scalar cwBeg = mag(blockPoints[v1] - blockPoints[v0]);
                const scalar cwEnd = mag(blockPoints[vn] - blockPoints[vn1]);

                Info<< "        j : "
                    << cwBeg*scaleTot.y() << " .. "
                    << cwEnd*scaleTot.y() << nl;
            }

            {
                const label nz = blocks[blocki].density().z();
                const label v1 = blocks[blocki].pointLabel(0, 0, 1);
                const label vn = blocks[blocki].pointLabel(0, 0, nz);
                const label vn1 = blocks[blocki].pointLabel(0, 0, nz-1);

                const scalar cwBeg = mag(blockPoints[v1] - blockPoints[v0]);
                const scalar cwEnd = mag(blockPoints[vn] - blockPoints[vn1]);

                Info<< "        k : "
                    << cwBeg*scaleTot.z() << " .. "
                    << cwEnd*scaleTot.z() << nl;
            }
            Info<< endl;
        }
    }

    inplacePointTransforms(points_);
}


void Foam::blockMesh::createCells() const
{
    const blockList& blocks = *this;
    const cellModel& hex = cellModel::ref(cellModel::HEX);

    if (verbose_)
    {
        Info<< "Creating cells" << endl;
    }

    cells_.resize(nCells_);

    label celli = 0;

    labelList cellPoints(8);  // Hex cells - 8 points

    forAll(blocks, blocki)
    {
        for (const hexCell& blockCell : blocks[blocki].cells())
        {
            forAll(cellPoints, cellPointi)
            {
                cellPoints[cellPointi] =
                    mergeList_
                    [
                        blockCell[cellPointi]
                      + blockOffsets_[blocki]
                    ];
            }

            // Construct collapsed cell and add to list
            cells_[celli].reset(hex, cellPoints, true);
            ++celli;
        }
    }
}


Foam::faceList Foam::blockMesh::createPatchFaces
(
    const polyPatch& patchTopologyFaces
) const
{
    const blockList& blocks = *this;

    labelList blockLabels = patchTopologyFaces.polyPatch::faceCells();

    label nFaces = 0;

    forAll(patchTopologyFaces, patchTopologyFaceLabel)
    {
        const label blocki = blockLabels[patchTopologyFaceLabel];

        faceList blockFaces = blocks[blocki].blockShape().faces();

        forAll(blockFaces, blockFaceLabel)
        {
            if
            (
                blockFaces[blockFaceLabel]
             == patchTopologyFaces[patchTopologyFaceLabel]
            )
            {
                nFaces +=
                    blocks[blocki].boundaryPatches()[blockFaceLabel].size();
            }
        }
    }


    faceList patchFaces(nFaces);
    face quad(4);
    label faceLabel = 0;

    forAll(patchTopologyFaces, patchTopologyFaceLabel)
    {
        const label blocki = blockLabels[patchTopologyFaceLabel];

        faceList blockFaces = blocks[blocki].blockShape().faces();

        forAll(blockFaces, blockFaceLabel)
        {
            if
            (
                blockFaces[blockFaceLabel]
             == patchTopologyFaces[patchTopologyFaceLabel]
            )
            {
                const List<FixedList<label, 4>>& blockPatchFaces =
                    blocks[blocki].boundaryPatches()[blockFaceLabel];

                forAll(blockPatchFaces, blockFaceLabel)
                {
                    // Lookup the face points
                    // and collapse duplicate point labels

                    quad[0] =
                        mergeList_
                        [
                            blockPatchFaces[blockFaceLabel][0]
                          + blockOffsets_[blocki]
                        ];

                    label nUnique = 1;

                    for
                    (
                        label facePointLabel = 1;
                        facePointLabel < 4;
                        facePointLabel++
                    )
                    {
                        quad[nUnique] =
                            mergeList_
                            [
                                blockPatchFaces[blockFaceLabel][facePointLabel]
                              + blockOffsets_[blocki]
                            ];

                        if (quad[nUnique] != quad[nUnique-1])
                        {
                            nUnique++;
                        }
                    }

                    if (quad[nUnique-1] == quad[0])
                    {
                        nUnique--;
                    }

                    if (nUnique == 4)
                    {
                        patchFaces[faceLabel++] = quad;
                    }
                    else if (nUnique == 3)
                    {
                        patchFaces[faceLabel++] = face
                        (
                            labelList::subList(quad, 3)
                        );
                    }
                    // else the face has collapsed to an edge or point
                }
            }
        }
    }

    patchFaces.resize(faceLabel);

    return patchFaces;
}


void Foam::blockMesh::createPatches() const
{
    const polyPatchList& topoPatches = topology().boundaryMesh();

    if (verbose_)
    {
        Info<< "Creating patches" << endl;
    }

    patches_.resize(topoPatches.size());

    forAll(topoPatches, patchi)
    {
        patches_[patchi] = createPatchFaces(topoPatches[patchi]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::blockMesh::topology() const
{
    if (!topologyPtr_)
    {
        FatalErrorInFunction
            << "topology not allocated"
            << exit(FatalError);
    }

    return *topologyPtr_;
}


Foam::refPtr<Foam::polyMesh>
Foam::blockMesh::topology(bool applyTransform) const
{
    const polyMesh& origTopo = topology();

    if (applyTransform && hasPointTransforms())
    {
        // Copy mesh components

        IOobject newIO(origTopo, IOobject::NO_READ, IOobject::NO_WRITE);
        newIO.registerObject(false);

        pointField newPoints(origTopo.points());
        inplacePointTransforms(newPoints);

        auto ttopoMesh = refPtr<polyMesh>::New
        (
            newIO,
            std::move(newPoints),
            faceList(origTopo.faces()),
            labelList(origTopo.faceOwner()),
            labelList(origTopo.faceNeighbour())
        );
        auto& topoMesh = ttopoMesh.ref();

        // Patches
        const polyBoundaryMesh& pbmOld = origTopo.boundaryMesh();
        const polyBoundaryMesh& pbmNew = topoMesh.boundaryMesh();

        polyPatchList newPatches(pbmOld.size());

        forAll(pbmOld, patchi)
        {
            newPatches.set(patchi, pbmOld[patchi].clone(pbmNew));
        }

        topoMesh.addPatches(newPatches);

        return ttopoMesh;
    }
    else
    {
        return origTopo;
    }
}


Foam::autoPtr<Foam::polyMesh>
Foam::blockMesh::mesh(const IOobject& io) const
{
    const blockMesh& blkMesh = *this;

    if (verbose_)
    {
        Info<< nl << "Creating polyMesh from blockMesh" << endl;
    }

    auto meshPtr = autoPtr<polyMesh>::New
    (
        io,
        pointField(blkMesh.points()),   // Use a copy of the block points
        blkMesh.cells(),
        blkMesh.patches(),
        blkMesh.patchNames(),
        blkMesh.patchDicts(),
        "defaultFaces",                 // Default patch name
        emptyPolyPatch::typeName        // Default patch type
    );


    // Set any cellZones
    const label nZones = blkMesh.numZonedBlocks();

    if (nZones)
    {
        polyMesh& pmesh = *meshPtr;

        if (verbose_)
        {
            Info<< "Adding cell zones" << endl;
        }

        // Map from zoneName to cellZone index
        HashTable<label> zoneMap(2*nZones);

        // Cells per zone
        List<DynamicList<label>> zoneCells(nZones);

        // Running cell counter
        label celli = 0;

        // Largest zone so far
        label freeZonei = 0;

        for (const block& b : blkMesh)
        {
            const word& zoneName = b.zoneName();
            const label nCellsInBlock = b.cells().size();

            if (zoneName.size())
            {
                const auto iter = zoneMap.cfind(zoneName);

                label zonei = freeZonei;

                if (iter.good())
                {
                    zonei = *iter;
                }
                else
                {
                    zoneMap.insert(zoneName, zonei);
                    ++freeZonei;

                    if (verbose_)
                    {
                        Info<< "    " << zonei << '\t' << zoneName << endl;
                    }
                }


                // Fill with cell ids

                zoneCells[zonei].reserve
                (
                    zoneCells[zonei].size() + nCellsInBlock
                );

                const label endOfFill = celli + nCellsInBlock;

                for (; celli < endOfFill; ++celli)
                {
                    zoneCells[zonei].append(celli);
                }
            }
            else
            {
                celli += nCellsInBlock;
            }
        }

        List<cellZone*> cz(zoneMap.size());
        forAllConstIters(zoneMap, iter)
        {
            const word& zoneName = iter.key();
            const label zonei = iter.val();

            cz[zonei] = new cellZone
            (
                zoneName,
                zoneCells[zonei].shrink(),
                zonei,
                pmesh.cellZones()
            );
        }

        pmesh.pointZones().clear();
        pmesh.faceZones().clear();
        pmesh.cellZones().clear();
        pmesh.addZones(List<pointZone*>(), List<faceZone*>(), cz);
    }


    // Merge patch pairs, cyclic must be done elsewhere
    // - requires libdynamicMesh

    return meshPtr;
}


// ************************************************************************* //
