/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "domainDecomposition.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "fvMesh.H"
#include "OSspecific.H"
#include "Map.H"
#include "DynamicList.H"
#include "fvFieldDecomposer.H"
#include "IOobjectList.H"
#include "PtrDynList.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "decompositionModel.H"
#include "hexRef8Data.H"

// For handling pointMeshes with additional patches
#include "pointMesh.H"
#include "meshPointPatch.H"
#include "processorPointPatch.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::domainDecomposition::mark
(
    const labelList& zoneElems,
    const label zoneI,
    labelList& elementToZone
)
{
    for (const label pointi : zoneElems)
    {
        if (elementToZone[pointi] == -1)
        {
            // First occurrence
            elementToZone[pointi] = zoneI;
        }
        else if (elementToZone[pointi] >= 0)
        {
            // Multiple zones
            elementToZone[pointi] = -2;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::domainDecomposition::domainDecomposition
(
    const IOobject& io,
    const fileName& decompDictFile
)
:
    fvMesh(io),
    facesInstancePointsPtr_
    (
        pointsInstance() != facesInstance()
      ? new pointIOField
        (
            IOobject
            (
                "points",
                facesInstance(),
                polyMesh::meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
      : nullptr
    ),
    decompDictFile_(decompDictFile),
    nProcs_
    (
        decompositionMethod::nDomains
        (
            decompositionModel::New
            (
                *this,
                decompDictFile
            )
        )
    ),
    distributed_(false),
    cellToProc_(nCells()),
    procPointAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
    procCellAddressing_(nProcs_),
    procPatchSize_(nProcs_),
    procPatchStartIndex_(nProcs_),
    procNeighbourProcessors_(nProcs_),
    procProcessorPatchSize_(nProcs_),
    procProcessorPatchStartIndex_(nProcs_),
    procProcessorPatchSubPatchIDs_(nProcs_),
    procProcessorPatchSubPatchStarts_(nProcs_)
{
    updateParameters(this->model());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::decompositionModel& Foam::domainDecomposition::model() const
{
    return decompositionModel::New(*this, decompDictFile_);
}


void Foam::domainDecomposition::updateParameters
(
    const dictionary& params
)
{
    params.readIfPresent("distributed", distributed_);
}


bool Foam::domainDecomposition::writeDecomposition(const bool decomposeSets)
{
    Info<< "\nConstructing processor meshes" << endl;

    // Mark point/faces/cells that are in zones.
    // -1   : not in zone
    // -2   : in multiple zones
    // >= 0 : in single given zone
    // This will give direct lookup of elements that are in a single zone
    // and we'll only have to revert back to searching through all zones
    // for the duplicate elements

    // Point zones
    labelList pointToZone(points().size(), -1);

    forAll(pointZones(), zonei)
    {
        mark(pointZones()[zonei], zonei, pointToZone);
    }

    // Face zones
    labelList faceToZone(faces().size(), -1);

    forAll(faceZones(), zonei)
    {
        mark(faceZones()[zonei], zonei, faceToZone);
    }

    // Cell zones
    labelList cellToZone(nCells(), -1);

    forAll(cellZones(), zonei)
    {
        mark(cellZones()[zonei], zonei, cellToZone);
    }


    PtrDynList<const cellSet> cellSets;
    PtrDynList<const faceSet> faceSets;
    PtrDynList<const pointSet> pointSets;
    if (decomposeSets)
    {
        // Read sets
        IOobjectList objects(*this, facesInstance(), "polyMesh/sets");
        for (const IOobject& io : objects.csorted<cellSet>())
        {
            cellSets.emplace_back(io);
        }
        for (const IOobject& io : objects.csorted<faceSet>())
        {
            faceSets.emplace_back(io);
        }
        for (const IOobject& io : objects.csorted<pointSet>())
        {
            pointSets.emplace_back(io);
        }
    }


    // Load refinement data (if any)
    hexRef8Data baseMeshData
    (
        IOobject
        (
            "dummy",
            facesInstance(),
            polyMesh::meshSubDir,
            *this,
            IOobjectOption::READ_IF_PRESENT,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        )
    );


    label maxProcCells = 0;
    label maxProcFaces = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;

    // Write out the meshes
    for (label proci = 0; proci < nProcs_; proci++)
    {
        // Create processor points
        const labelList& curPointLabels = procPointAddressing_[proci];

        const pointField& meshPoints = points();

        labelList pointLookup(nPoints(), -1);

        pointField procPoints(curPointLabels.size());

        forAll(curPointLabels, pointi)
        {
            procPoints[pointi] = meshPoints[curPointLabels[pointi]];

            pointLookup[curPointLabels[pointi]] = pointi;
        }

        // Create processor faces
        const labelList& curFaceLabels = procFaceAddressing_[proci];

        const faceList& meshFaces = faces();

        labelList faceLookup(nFaces(), -1);

        faceList procFaces(curFaceLabels.size());

        forAll(curFaceLabels, facei)
        {
            // Mark the original face as used
            // Remember to decrement the index by one (turning index)
            label curF = mag(curFaceLabels[facei]) - 1;

            faceLookup[curF] = facei;

            // get the original face
            labelList origFaceLabels;

            if (curFaceLabels[facei] >= 0)
            {
                // face not turned
                origFaceLabels = meshFaces[curF];
            }
            else
            {
                origFaceLabels = meshFaces[curF].reverseFace();
            }

            // translate face labels into local point list
            face& procFaceLabels = procFaces[facei];

            procFaceLabels.setSize(origFaceLabels.size());

            forAll(origFaceLabels, pointi)
            {
                procFaceLabels[pointi] = pointLookup[origFaceLabels[pointi]];
            }
        }

        // Create processor cells
        const labelList& curCellLabels = procCellAddressing_[proci];

        const cellList& meshCells = cells();

        cellList procCells(curCellLabels.size());

        forAll(curCellLabels, celli)
        {
            const labelList& origCellLabels = meshCells[curCellLabels[celli]];

            cell& curCell = procCells[celli];

            curCell.setSize(origCellLabels.size());

            forAll(origCellLabels, cellFacei)
            {
                curCell[cellFacei] = faceLookup[origCellLabels[cellFacei]];
            }
        }

        // Create processor mesh without a boundary

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(proci)),
            false,  // No function objects
            false   // No extra controlDict libs
        );
        processorDb.setTime(time());

        // create the mesh. Two situations:
        // - points and faces come from the same time ('instance'). The mesh
        //   will get constructed in the same instance.
        // - points come from a different time (moving mesh cases).
        //   It will read the points belonging to the faces instance and
        //   construct the procMesh with it which then gets handled as above.
        //   (so with 'old' geometry).
        //   Only at writing time will it additionally write the current
        //   points.

        autoPtr<polyMesh> procMeshPtr;

        if (facesInstancePointsPtr_)
        {
            // Construct mesh from facesInstance.
            pointField facesInstancePoints
            (
                facesInstancePointsPtr_(),
                curPointLabels
            );

            procMeshPtr = autoPtr<polyMesh>::New
            (
                IOobject
                (
                    this->polyMesh::name(), // region of undecomposed mesh
                    facesInstance(),
                    processorDb,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                std::move(facesInstancePoints),
                std::move(procFaces),
                std::move(procCells)
            );
        }
        else
        {
            procMeshPtr = autoPtr<polyMesh>::New
            (
                IOobject
                (
                    this->polyMesh::name(), // region of undecomposed mesh
                    facesInstance(),
                    processorDb,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                std::move(procPoints),
                std::move(procFaces),
                std::move(procCells)
            );
        }
        polyMesh& procMesh = procMeshPtr();


        // Create processor boundary patches
        const labelList& curPatchSizes = procPatchSize_[proci];

        const labelList& curPatchStarts = procPatchStartIndex_[proci];

        const labelList& curNeighbourProcessors =
            procNeighbourProcessors_[proci];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[proci];

        const labelList& curProcessorPatchStarts =
            procProcessorPatchStartIndex_[proci];

        const labelListList& curSubPatchIDs =
            procProcessorPatchSubPatchIDs_[proci];

        const labelListList& curSubStarts =
            procProcessorPatchSubPatchStarts_[proci];

        const polyPatchList& meshPatches = boundaryMesh();


        // Count the number of inter-proc patches
        label nInterProcPatches = 0;
        forAll(curSubPatchIDs, procPatchi)
        {
            nInterProcPatches += curSubPatchIDs[procPatchi].size();
        }

        polyPatchList procPatches
        (
            curPatchSizes.size() + nInterProcPatches
        );

        label nPatches = 0;

        forAll(curPatchSizes, patchi)
        {
            // Get the face labels consistent with the field mapping
            // (reuse the patch field mappers)
            const polyPatch& meshPatch = meshPatches[patchi];

            fvFieldDecomposer::patchFieldDecomposer patchMapper
            (
                SubList<label>
                (
                    curFaceLabels,
                    curPatchSizes[patchi],
                    curPatchStarts[patchi]
                ),
                meshPatch.start()
            );

            // Map existing patches
            procPatches.set
            (
                nPatches,
                meshPatch.clone
                (
                    procMesh.boundaryMesh(),
                    nPatches,
                    patchMapper.directAddressing(),
                    curPatchStarts[patchi]
                )
            );

            nPatches++;
        }

        forAll(curProcessorPatchSizes, procPatchi)
        {
            const labelList& subPatchID = curSubPatchIDs[procPatchi];
            const labelList& subStarts = curSubStarts[procPatchi];

            label curStart = curProcessorPatchStarts[procPatchi];

            forAll(subPatchID, i)
            {
                label size =
                (
                    i < subPatchID.size()-1
                  ? subStarts[i+1] - subStarts[i]
                  : curProcessorPatchSizes[procPatchi] - subStarts[i]
                );

                if (subPatchID[i] == -1)
                {
                    // From internal faces
                    procPatches.set
                    (
                        nPatches,
                        new processorPolyPatch
                        (
                            size,
                            curStart,
                            nPatches,
                            procMesh.boundaryMesh(),
                            proci,
                            curNeighbourProcessors[procPatchi]
                        )
                    );
                }
                else
                {
                    const coupledPolyPatch& pcPatch
                        = refCast<const coupledPolyPatch>
                          (
                              boundaryMesh()[subPatchID[i]]
                          );

                    procPatches.set
                    (
                        nPatches,
                        new processorCyclicPolyPatch
                        (
                            size,
                            curStart,
                            nPatches,
                            procMesh.boundaryMesh(),
                            proci,
                            curNeighbourProcessors[procPatchi],
                            pcPatch.name(),
                            pcPatch.transform()
                        )
                    );
                }

                curStart += size;
                ++nPatches;
            }
        }

        // Add boundary patches
        procMesh.addPatches(procPatches);

        // Create and add zones

        // Point zones
        {
            const pointZoneMesh& pz = pointZones();

            // Go through all the zoned points and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zonePoints(pz.size());

            // Estimate size
            forAll(zonePoints, zonei)
            {
                zonePoints[zonei].setCapacity(pz[zonei].size()/nProcs_);
            }

            // Use the pointToZone map to find out the single zone (if any),
            // use slow search only for shared points.
            forAll(curPointLabels, pointi)
            {
                label curPoint = curPointLabels[pointi];

                label zonei = pointToZone[curPoint];

                if (zonei >= 0)
                {
                    // Single zone.
                    zonePoints[zonei].append(pointi);
                }
                else if (zonei == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(pz, zonei)
                    {
                        label index = pz[zonei].whichPoint(curPoint);

                        if (index != -1)
                        {
                            zonePoints[zonei].append(pointi);
                        }
                    }
                }
            }

            procMesh.pointZones().clearAddressing();
            procMesh.pointZones().setSize(zonePoints.size());
            forAll(zonePoints, zonei)
            {
                procMesh.pointZones().set
                (
                    zonei,
                    pz[zonei].clone
                    (
                        procMesh.pointZones(),
                        zonei,
                        zonePoints[zonei].shrink()
                    )
                );
            }

            if (pz.size())
            {
                // Force writing on all processors
                procMesh.pointZones().writeOpt(IOobject::AUTO_WRITE);
            }
        }

        // Face zones
        {
            const faceZoneMesh& fz = faceZones();

            // Go through all the zoned face and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zoneFaces(fz.size());
            List<DynamicList<bool>> zoneFaceFlips(fz.size());

            // Estimate size
            forAll(zoneFaces, zonei)
            {
                label procSize = fz[zonei].size() / nProcs_;

                zoneFaces[zonei].setCapacity(procSize);
                zoneFaceFlips[zonei].setCapacity(procSize);
            }

            // Go through all the zoned faces and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            forAll(curFaceLabels, facei)
            {
                // Remember to decrement the index by one (turning index)
                //
                label curF = mag(curFaceLabels[facei]) - 1;

                label zonei = faceToZone[curF];

                if (zonei >= 0)
                {
                    // Single zone. Add the face
                    zoneFaces[zonei].append(facei);

                    label index = fz[zonei].whichFace(curF);

                    bool flip = fz[zonei].flipMap()[index];

                    if (curFaceLabels[facei] < 0)
                    {
                        flip = !flip;
                    }

                    zoneFaceFlips[zonei].append(flip);
                }
                else if (zonei == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(fz, zonei)
                    {
                        label index = fz[zonei].whichFace(curF);

                        if (index != -1)
                        {
                            zoneFaces[zonei].append(facei);

                            bool flip = fz[zonei].flipMap()[index];

                            if (curFaceLabels[facei] < 0)
                            {
                                flip = !flip;
                            }

                            zoneFaceFlips[zonei].append(flip);
                        }
                    }
                }
            }

            procMesh.faceZones().clearAddressing();
            procMesh.faceZones().setSize(zoneFaces.size());
            forAll(zoneFaces, zonei)
            {
                procMesh.faceZones().set
                (
                    zonei,
                    fz[zonei].clone
                    (
                        zoneFaces[zonei].shrink(),          // addressing
                        zoneFaceFlips[zonei].shrink(),      // flipmap
                        zonei,
                        procMesh.faceZones()
                    )
                );
            }

            if (fz.size())
            {
                // Force writing on all processors
                procMesh.faceZones().writeOpt(IOobject::AUTO_WRITE);
            }
        }

        // Cell zones
        {
            const cellZoneMesh& cz = cellZones();

            // Go through all the zoned cells and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zoneCells(cz.size());

            // Estimate size
            forAll(zoneCells, zonei)
            {
                zoneCells[zonei].setCapacity(cz[zonei].size()/nProcs_);
            }

            forAll(curCellLabels, celli)
            {
                label curCelli = curCellLabels[celli];

                label zonei = cellToZone[curCelli];

                if (zonei >= 0)
                {
                    // Single zone.
                    zoneCells[zonei].append(celli);
                }
                else if (zonei == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(cz, zonei)
                    {
                        label index = cz[zonei].whichCell(curCelli);

                        if (index != -1)
                        {
                            zoneCells[zonei].append(celli);
                        }
                    }
                }
            }

            procMesh.cellZones().clearAddressing();
            procMesh.cellZones().setSize(zoneCells.size());
            forAll(zoneCells, zonei)
            {
                procMesh.cellZones().set
                (
                    zonei,
                    cz[zonei].clone
                    (
                        zoneCells[zonei].shrink(),
                        zonei,
                        procMesh.cellZones()
                    )
                );
            }

            if (cz.size())
            {
                // Force writing on all processors
                procMesh.cellZones().writeOpt(IOobject::AUTO_WRITE);
            }
        }

        // More precision (for points data)
        IOstream::minPrecision(10);

        procMesh.write();

        // Add pointMesh if it was available
        const auto* pMeshPtr =
            thisDb().cfindObject<pointMesh>(pointMesh::typeName);
        if (pMeshPtr)
        {
            const auto& pMesh = *pMeshPtr;
            const auto& pMeshBoundary = pMesh.boundary();


            // 1. Generate pointBoundaryMesh from polyBoundaryMesh (so ignoring
            //    any additional patches
            const auto& procPointMesh = pointMesh::New(procMesh);

            pointBoundaryMesh& procBoundary =
                const_cast<pointBoundaryMesh&>(procPointMesh.boundary());

            // Keep track if it differs from the polyBoundaryMesh since then
            // we need to write the boundary file.
            bool differsFromPoly = false;

            // 2. Explicitly add subsetted meshPointPatches
            forAll(pMeshBoundary, patchi)
            {
                const auto* mppPtr = isA<meshPointPatch>(pMeshBoundary[patchi]);
                if (mppPtr && (procBoundary.findPatchID(mppPtr->name()) == -1))
                {
                    const auto& mpp = *mppPtr;

                    DynamicList<label> procMeshPoints(mpp.size());
                    DynamicField<vector> procNormals(mpp.size());
                    forAll(mpp.meshPoints(), i)
                    {
                        const label pointi = mpp.meshPoints()[i];
                        const label procPointi = pointLookup[pointi];
                        if (procPointi != -1)
                        {
                            procMeshPoints.append(procPointi);
                            procNormals.append(mpp.pointNormals()[i]);
                        }
                    }

                    procBoundary.push_back
                    (
                        new meshPointPatch
                        (
                            mpp.name(),
                            procMeshPoints,
                            procNormals,
                            procBoundary.size(),
                            procBoundary,
                            meshPointPatch::typeName
                        )
                    );
                    differsFromPoly = true;
                }
            }

            // 3. Shuffle new patches before any processor patches
            labelList oldToNew(procBoundary.size());
            label newPatchi = 0;
            forAll(procBoundary, patchi)
            {
                if (!isA<processorPointPatch>(procBoundary[patchi]))
                {
                    oldToNew[patchi] = newPatchi;

                    if (newPatchi != patchi)
                    {
                        differsFromPoly = true;
                    }

                    newPatchi++;
                }
            }

            // decomposed-to-undecomposed patch numbering
            labelList boundaryProcAddressing(identity(newPatchi));
            boundaryProcAddressing.setSize(procBoundary.size(), -1);

            forAll(procBoundary, patchi)
            {
                if (isA<processorPointPatch>(procBoundary[patchi]))
                {
                    oldToNew[patchi] = newPatchi++;
                }
            }
            procBoundary.reorder(oldToNew, true);


            if (differsFromPoly)
            {
                // Write pointMesh/boundary
                procBoundary.write();

                // Write pointMesh/boundaryProcAddressing
                IOobject ioAddr
                (
                    "boundaryProcAddressing",
                    procMesh.facesInstance(),
                    polyMesh::meshSubDir/pointMesh::meshSubDir,
                    procPointMesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                );
                IOListRef<label>(ioAddr, boundaryProcAddressing).write();
            }
        }

        // Write points if pointsInstance differing from facesInstance
        if (facesInstancePointsPtr_)
        {
            pointIOField pointsInstancePoints
            (
                IOobject
                (
                    "points",
                    pointsInstance(),
                    polyMesh::meshSubDir,
                    procMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                std::move(procPoints)
            );
            pointsInstancePoints.write();
        }


        // Decompose any sets
        if (decomposeSets)
        {
            forAll(cellSets, i)
            {
                const cellSet& cs = cellSets[i];
                cellSet set(procMesh, cs.name(), cs.size()/nProcs_);
                forAll(curCellLabels, i)
                {
                    if (cs.found(curCellLabels[i]))
                    {
                        set.insert(i);
                    }
                }
                set.write();
            }
            forAll(faceSets, i)
            {
                const faceSet& cs = faceSets[i];
                faceSet set(procMesh, cs.name(), cs.size()/nProcs_);
                forAll(curFaceLabels, i)
                {
                    if (cs.found(mag(curFaceLabels[i])-1))
                    {
                        set.insert(i);
                    }
                }
                set.write();
            }
            forAll(pointSets, i)
            {
                const pointSet& cs = pointSets[i];
                pointSet set(procMesh, cs.name(), cs.size()/nProcs_);
                forAll(curPointLabels, i)
                {
                    if (cs.found(curPointLabels[i]))
                    {
                        set.insert(i);
                    }
                }
                set.write();
            }
        }


        // Optional hexRef8 data
        hexRef8Data
        (
            IOobject
            (
                "dummy",
                facesInstance(),
                polyMesh::meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            baseMeshData,
            procCellAddressing_[proci],
            procPointAddressing_[proci]
        ).write();


        // Statistics
        Info<< nl << "Processor " << proci;

        if (procMesh.nCells())
        {
            Info<< nl << "    ";
        }
        else
        {
            Info<< ": ";
        }

        Info<< "Number of cells = " << procMesh.nCells() << nl;

        if (procMesh.nCells())
        {
            Info<< "    Number of points = " << procMesh.nPoints() << nl;
        }

        maxProcCells = max(maxProcCells, procMesh.nCells());

        label nBoundaryFaces = 0;
        label nProcPatches = 0;
        label nProcFaces = 0;

        for (const polyPatch& pp : procMesh.boundaryMesh())
        {
            const auto* cpp = isA<processorPolyPatch>(pp);

            if (cpp)
            {
                const auto& procPatch = *cpp;

                Info<< "    Number of faces shared with processor "
                    << procPatch.neighbProcNo() << " = "
                    << procPatch.size() << nl;

                nProcFaces += procPatch.size();
                ++nProcPatches;
            }
            else
            {
                nBoundaryFaces += pp.size();
            }
        }

        if (procMesh.nCells() && (nBoundaryFaces || nProcFaces))
        {
            Info<< "    Number of processor patches = " << nProcPatches << nl
                << "    Number of processor faces = " << nProcFaces << nl
                << "    Number of boundary faces = " << nBoundaryFaces << nl;
        }

        totProcFaces += nProcFaces;
        totProcPatches += nProcPatches;
        maxProcFaces = max(maxProcFaces, nProcFaces);
        maxProcPatches = max(maxProcPatches, nProcPatches);

        // Write the addressing information

        IOobject ioAddr
        (
            "procAddressing",
            procMesh.facesInstance(),
            polyMesh::meshSubDir,
            procMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        // pointProcAddressing
        ioAddr.rename("pointProcAddressing");
        IOListRef<label>(ioAddr, procPointAddressing_[proci]).write();

        // faceProcAddressing
        ioAddr.rename("faceProcAddressing");
        IOListRef<label>(ioAddr, procFaceAddressing_[proci]).write();

        // cellProcAddressing
        ioAddr.rename("cellProcAddressing");
        IOListRef<label>(ioAddr, procCellAddressing_[proci]).write();

        // Write patch map for backwards compatibility.
        // (= identity map for original patches, -1 for processor patches)
        label nMeshPatches = curPatchSizes.size();
        labelList procBoundaryAddr(identity(nMeshPatches));
        procBoundaryAddr.resize(nMeshPatches+nProcPatches, -1);

        // boundaryProcAddressing
        ioAddr.rename("boundaryProcAddressing");
        IOListRef<label>(ioAddr, procBoundaryAddr).write();
    }


    // Summary stats
    Info<< nl
        << "Number of processor faces = " << (totProcFaces/2) << nl
        << "Max number of cells = " << maxProcCells;

    if (maxProcCells != nCells())
    {
        scalar avgValue = scalar(nCells())/nProcs_;

        Info<< " (" << 100.0*(maxProcCells-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl;

    Info<< "Max number of processor patches = " << maxProcPatches;
    if (totProcPatches)
    {
        scalar avgValue = scalar(totProcPatches)/nProcs_;

        Info<< " (" << 100.0*(maxProcPatches-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl;

    Info<< "Max number of faces between processors = " << maxProcFaces;
    if (totProcFaces)
    {
        scalar avgValue = scalar(totProcFaces)/nProcs_;

        Info<< " (" << 100.0*(maxProcFaces-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl << endl;

    return true;
}


// ************************************************************************* //
