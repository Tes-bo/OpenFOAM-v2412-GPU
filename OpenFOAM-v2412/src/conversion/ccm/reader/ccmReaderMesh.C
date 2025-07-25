/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "ccmReader.H"

#include "emptyPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wallPolyPatch.H"
#include "Fstream.H"
#include "IOdictionary.H"

#include "ccmBoundaryInfo.H"
#include "indirectPrimitivePatch.H"
#include "SortableList.H"
#include "mergePoints.H"
#include "bitSet.H"
#include "ListOps.H"

#include "ccmInternal.H" // include last to avoid any strange interactions

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::ccm::reader::patchStartList(label initial) const
{
    labelList startLst(patchSizes_.size(), Zero);

    label patchFaceI = initial;
    forAll(patchSizes_, patchI)
    {
        startLst[patchI] = patchFaceI;
        patchFaceI += patchSizes_[patchI];
    }

    return startLst;
}


void Foam::ccm::reader::printSizes() const
{
    Info<<"nPoints:" << nPoints_
        << "  nCells:" << nCells_
        << "  nFaces:" << nFaces_
        << "  nInternalFaces:" << nInternalFaces_
        << endl;
}


// Determine if the geometry looks good.
// We'll insist that everything be on the first ("default") state,
// and is on the same file, and has a problem description
//
// The detection reads the problem description as well
//
bool Foam::ccm::reader::detectGeometry()
{
    // Call once
    if (geometryStatus_ != UNKNOWN)
    {
        return (geometryStatus_ == OKAY || geometryStatus_ == READ);
    }

    // Explicitly restricted to 'default' state node
    ccmID stateNode;

    // Get associated problem description
    ccmID probNode;

    // Only check the first processor
    ccmID processorNode;
    int procI = 0;

    // Geometry needs vertices and topology
    ccmID verticesNode, topoNode;

    if
    (
        CCMIOGetState
        (
            nullptr,
            (globalState_->root),
            "default",
            &probNode,
            &stateNode
        )
     == kCCMIONoErr

     && CCMIONextEntity
        (
            nullptr,
            stateNode,
            kCCMIOProcessor,
            &procI,
            &processorNode
        )
     == kCCMIONoErr

     && CCMIOReadProcessor
        (
            nullptr,
            processorNode,
            &verticesNode,
            &topoNode,
            nullptr,        // Ignore initialField
            nullptr         // Ignore solutionNode
        )
     == kCCMIONoErr

     && CCMIOIsFromSameFile((globalState_->root), verticesNode)
     && CCMIOIsFromSameFile((globalState_->root), topoNode)
     && CCMIOIsValidEntity(probNode)
    )
    {
        readProblemDescription(probNode);
        geometryStatus_ = OKAY;
    }
    else
    {
        // Missing/incomplete geometry node and/or problem node
        geometryStatus_ = BAD;
    }

    return (geometryStatus_ == OKAY || geometryStatus_ == READ);
}


// read the geometry without any sorting
void Foam::ccm::reader::readMeshTopology
(
    const scalar scaleFactor
)
{
    ccmID verticesNode, topoNode;

    // Use first ("default") state node
    ccmID stateNode;
    int stateI = 0;

    // Use the first processor to find the mesh nodes
    ccmID processorNode;
    int procI = 0;

    if
    (
        CCMIONextEntity
        (
            nullptr,
            (globalState_->root),
            kCCMIOState,
            &stateI,
            &stateNode
        )
     == kCCMIONoErr

     && CCMIONextEntity
        (
            nullptr,
            stateNode,
            kCCMIOProcessor,
            &procI,
            &processorNode
        )
     == kCCMIONoErr

     && CCMIOReadProcessor
        (
            nullptr,
            processorNode,
            &verticesNode,
            &topoNode,
            nullptr,        // Ignore initialField
            nullptr         // Ignore solutionNode
        )
     == kCCMIONoErr
    )
    {
        labelList origPointId = readVertices(verticesNode, scaleFactor);
        readCells(topoNode);
        readMonitoring(topoNode);

        // Renumber vertex labels (ccm -> Foam)
        {
            label maxId = max(origPointId);
            labelList mapToFoam(invert(maxId+1, origPointId));

            for (face& f : faces_)
            {
                inplaceRenumber(mapToFoam, f);
            }
        }
        origPointId.clear();

        // Renumber owners/neighbours cell labels (ccm -> Foam)
        {
            label maxId = max(origCellId_);
            labelList mapToFoam(invert(maxId+1, origCellId_));

            inplaceRenumber(mapToFoam, faceOwner_);
            inplaceRenumber(mapToFoam, faceNeighbour_);
        }

        // Juggle solids into fluid as required
        juggleSolids();

        // Report sizes
        printSizes();

        // Remove unwanted fluid/porous/solid types
        removeUnwanted();

        // Collapse interfaces between domains (eg, fluid|porosity)
        cleanupInterfaces();

        // Use point merge to join conformal interfaces (STARCCM)
        mergeInplaceInterfaces();
    }
}


// readVertices:
// 1) read the vertex data
// 2) read the map (which maps the index into the data array with the Id number)
//
// returns the original point Id
Foam::labelList Foam::ccm::reader::readVertices
(
    const ccmID& verticesNode,
    const scalar scaleFactor
)
{
    int dims = 1;
    float scale = 1.0;

    CCMIOSize size = 0;
    ccmID mapId;

#ifdef DEBUG_CCMIOREAD
    Info<< "readVertices()" << endl;
#endif

    // Determine dimensions and mapId
    CCMIOEntitySize
    (
        &(globalState_->error),
        verticesNode,
        &size,
        nullptr
    );

    CCMIOReadVerticesd
    (
        &(globalState_->error),
        verticesNode,
        &dims,
        &scale,
        &mapId,
        nullptr,
        kCCMIOStart,
        kCCMIOEnd
    );
    assertNoError("problem finding 'Vertices' node");

    if (dims != 3)
    {
        FatalErrorInFunction
            << "can only handle 3-dimensional vertices"
            << exit(FatalError);
    }

    nPoints_ = size;
    labelList origPointId(nPoints_);

    readMap
    (
        mapId,
        origPointId
    );

    // Temporary storage for the points - reading piecemeal is much slower
    List<scalar> vrts(3*nPoints_);

    // The ccm data is double precision too
    CCMIOReadVerticesd
    (
        &(globalState_->error),
        verticesNode,
        nullptr,
        nullptr,
        nullptr,
        vrts.data(),
        kCCMIOStart,
        kCCMIOEnd
    );
    assertNoError("problem reading 'Vertices' node");

    // Convert to foam Points
    points_.setSize(nPoints_);

    scalar effectiveScale = scale * scaleFactor;
    forAll(points_, i)
    {
        points_[i].x() = effectiveScale * vrts[i*3];
        points_[i].y() = effectiveScale * vrts[i*3+1];
        points_[i].z() = effectiveScale * vrts[i*3+2];
    }

    vrts.clear();

#ifdef DEBUG_CCMIOREAD
    Info<< "readVertices: " << nPoints_ << endl;
#endif

    return origPointId;
}


// readCells:
// - read faces, faceOwner and faceNeighbour
// - finally read interfaces
void Foam::ccm::reader::readCells
(
    const ccmID& topoNode
)
{
    CCMIOSize size = 0;
    ccmID cellsNode, mapId;
    ccmID nodeId;

#ifdef DEBUG_CCMIOREAD
    Info<< "readCells()" << endl;
#endif

    // Determine dimensions and mapId information for 'Cells'
    CCMIOGetEntity
    (
        &(globalState_->error),
        topoNode,
        kCCMIOCells,
        0,
        &cellsNode
    );

    CCMIOEntitySize
    (
        &(globalState_->error),
        cellsNode,
        &size,
        nullptr
    );
    assertNoError("cannot get 'Cells' node");

    nCells_ = size;

#ifdef DEBUG_CCMIOREAD
    Info<< "readCells: " << nCells_ << endl;
#endif

    // Store cell ids so that we know which cells are which
    origCellId_.setSize(nCells_);
    cellTableId_.setSize(nCells_);

    CCMIOReadCells
    (
        &(globalState_->error),
        cellsNode,
        &mapId,
        cellTableId_.data(),
        kCCMIOStart,
        kCCMIOEnd
    );
    assertNoError("Error reading 'Cells' node");

    readMap
    (
        mapId,
        origCellId_
    );

    // Determine dimensions and mapId information for 'InternalFaces'
    CCMIOGetEntity
    (
        &(globalState_->error),
        topoNode,
        kCCMIOInternalFaces,
        0,
        &nodeId
    );
    CCMIOEntitySize
    (
        &(globalState_->error),
        nodeId,
        &size,
        nullptr
    );
    nInternalFaces_ = size;

    nFaces_ = nInternalFaces_;

    // First pass:
    //
    // Determine patch sizes before reading internal faces
    // also determine the original boundary regions

    label nPatches = 0;
    DynamicList<ccmBoundaryInfo> bndInfo;

    // Number of children in the parent node is more than number of patches,
    // but is a good start for allocation
    if
    (
        CCMIOGetNumberOfChildren
        (
            nullptr,
            topoNode.node,
            &nPatches
        ) == kCCMIONoErr
    )
    {
        bndInfo.setCapacity(nPatches);
    }

    for
    (
        int index = 0;
        CCMIONextEntity
        (
            nullptr,
            topoNode,
            kCCMIOBoundaryFaces,
            &index,
            &nodeId
        ) == kCCMIONoErr
     && CCMIOEntitySize
        (
            &(globalState_->error),
            nodeId,
            &size,
            nullptr
        ) == kCCMIONoErr;
        /* nop */
    )
    {
        ccmBoundaryInfo info;
        info.size = size;
        nFaces_   += size;

        CCMIOGetEntityIndex
        (
            &(globalState_->error),
            nodeId,
            &(info.ccmIndex)
        );

        // Name directly from the node (eg, STARCCM)
        info.setPatchName(ccmReadOptstr("Label", nodeId));

        // Lookup the name, type from boundary region info:
        auto dictIter = boundaryRegion_.find(info.ccmIndex);
        if (dictIter.good())
        {
            dictionary& dict = dictIter.val();

            const word patchName(dict.get<word>("Label"));
            const word patchType(dict.get<word>("BoundaryType"));

            if (!patchName.empty())
            {
                info.patchName = patchName;
            }

            if (!patchType.empty())
            {
                info.patchType = patchType;
            }

            // Optional, but potentially useful information:
            dict.add("BoundaryIndex", info.ccmIndex);
            dict.add("size", info.size);
        }

        bndInfo.append(info);
    }

    // Redimension lists according to the overall sizes
    faces_.setSize(nFaces_);
    faceOwner_.setSize(nFaces_);
    faceNeighbour_.setSize(nFaces_);
    origFaceId_.setSize(nFaces_);


    // May be too large, but is a good place start size
    patchSizes_.setSize(bndInfo.size());
    origBndId_.setSize(bndInfo.size());
    patchSizes_ = 0;
    origBndId_  = -1;
    nPatches = 0;

    for (ccmBoundaryInfo& info : bndInfo)
    {
        if (info.patchId != -1)
        {
            // Already inserted - eg, as interface
            origBndId_[info.patchId] = info.ccmIndex;
        }
        else
        {
            info.patchId = nPatches++;
            origBndId_[info.patchId] = info.ccmIndex;
        }

        patchSizes_[info.patchId] += info.size;
    }

    // Shrink to sizes actually used
    patchSizes_.setSize(nPatches);
    origBndId_.setSize(nPatches);

    // Boundary info indices flattened and sorted by patchId

    IndirectList<ccmBoundaryInfo> ccmLookupOrder(bndInfo, labelList());
    {
        DynamicList<label> addr(bndInfo.size());
        for (int patchI = 0; patchI < nPatches; ++patchI)
        {
            forAll(bndInfo, infoI)
            {
                if (bndInfo[infoI].patchId == patchI)
                {
                    addr.append(infoI);
                }
            }
        }

        ccmLookupOrder.addressing() = std::move(addr);
    }


    //
    // Now we are ready to do the reading
    //
    CCMIOGetEntity
    (
        &(globalState_->error),
        topoNode,
        kCCMIOInternalFaces,
        0,
        &nodeId
    );

    // get allocation sizes
    CCMIOReadFaces
    (
        &(globalState_->error),
        nodeId,
        kCCMIOInternalFaces,
        &mapId,
        &size,
        nullptr,
        kCCMIOStart,
        kCCMIOEnd
    );

    List<label> mapData(nInternalFaces_);
    List<label> faceCells(2*nInternalFaces_);
    List<label> ccmFaces(size);

    CCMIOReadFaces
    (
        &(globalState_->error),
        nodeId,
        kCCMIOInternalFaces,
        nullptr,
        nullptr,
        ccmFaces.data(),
        kCCMIOStart,
        kCCMIOEnd
    );

    CCMIOReadFaceCells
    (
        &(globalState_->error),
        nodeId,
        kCCMIOInternalFaces,
        faceCells.data(),
        kCCMIOStart,
        kCCMIOEnd
    );
    assertNoError("Error reading internal faces");

    readMap
    (
        mapId,
        mapData
    );

    // Copy into Foam list
    // ccmFaces are organized as [nVert vrt1 .. vrtN]
    unsigned int pos = 0;
    for (label faceI = 0; faceI < nInternalFaces_; ++faceI)
    {
        origFaceId_[faceI]    = mapData[faceI];
        faceOwner_[faceI]     = faceCells[2*faceI];
        faceNeighbour_[faceI] = faceCells[2*faceI+1];
        face& f = faces_[faceI];

        f.setSize(ccmFaces[pos++]);
        forAll(f, fp)
        {
            f[fp] = ccmFaces[pos++];
        }
    }

    // Read the boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~
    label patchFaceI = nInternalFaces_;

    for (const ccmBoundaryInfo& info : ccmLookupOrder)
    {
        const unsigned int patchSize = info.size;

        if
        (
            CCMIOGetEntity
            (
                &(globalState_->error),
                topoNode,
                kCCMIOBoundaryFaces,
                info.ccmIndex,
                &nodeId
            ) == kCCMIONoErr
        )
        {
            CCMIOReadFaces
            (
                &(globalState_->error),
                nodeId,
                kCCMIOBoundaryFaces,
                &mapId,
                &size,    // size needed for the faces
                nullptr,
                kCCMIOStart,
                kCCMIOEnd
            );

            mapData.setSize(patchSize);
            faceCells.setSize(patchSize);
            ccmFaces.setSize(size);

            readMap
            (
                mapId,
                mapData
            );

            CCMIOReadFaces
            (
                &(globalState_->error),
                nodeId,
                kCCMIOBoundaryFaces,
                nullptr,
                nullptr,
                ccmFaces.data(),
                kCCMIOStart,
                kCCMIOEnd
            );
            CCMIOReadFaceCells
            (
                &(globalState_->error),
                nodeId,
                kCCMIOBoundaryFaces,
                faceCells.data(),
                kCCMIOStart,
                kCCMIOEnd
            );
            assertNoError
            (
                "Error reading boundary face cells - index "
              + ::Foam::name(info.ccmIndex)
            );

            // Copy into Foam list
            // ccmFaces are organized as [nVert vrt1 .. vrtN]
            unsigned int pos = 0;
            for (unsigned int i = 0; i < patchSize; ++i, ++patchFaceI)
            {
                origFaceId_[patchFaceI]    = mapData[i];
                faceOwner_[patchFaceI]     = faceCells[i];
                faceNeighbour_[patchFaceI] = -1;

                face& f = faces_[patchFaceI];

                f.setSize(ccmFaces[pos++]);
                forAll(f, fp)
                {
                    f[fp] = ccmFaces[pos++];
                }
            }
        }
        else
        {
            assertNoError
            (
                "Error reading boundary faces - index "
              + ::Foam::name(info.ccmIndex)
            );
        }
    }

    readInterfaces(cellsNode);
}


// Read any interfaces
//   1) interfaces between domains (fluid/solid, fluid/porosity)
//   2) PROSTAR baffles
//
void Foam::ccm::reader::readInterfaces
(
    const ccmID& cellsNode
)
{
#ifdef DEBUG_CCMIOREAD
    Info<< "readInterfaces()" << endl;
#endif

    label nBaffleInterface = 0, nInterfaceTotal = 0;
    CCMIOIndex size = 0, dims = 0;

    bafInterfaces_.clear();
    domInterfaces_.clear();

    ccmID interfaceNode;

    if
    (
        CCMIOGetEntity
        (
            nullptr,
            cellsNode,
            kCCMIOInterfaces,
            0,
            &interfaceNode
        )
     == kCCMIONoErr

     && CCMIOGetOptInfo
        (
            nullptr,
            interfaceNode,
            "FaceIds",
            nullptr,
            &size,
            &dims,
            nullptr
        )
     == kCCMIONoErr

     && size > 0 && dims == 2
    )
    {
        nInterfaceTotal = size;
    }
    else
    {
        return;
    }

    // Get ProstarBaffles
    // formatted as [ prostarId faceId1 faceId2 celltableId ]
    if
    (
        CCMIOGetOptInfo
        (
            nullptr,
            interfaceNode,
            "ProstarBaffles",
            nullptr,
            &size,
            nullptr,
            nullptr
        )
     == kCCMIONoErr

     && size > 0
    )
    {
        // Be paranoid - force an integral value of 4
        nBaffleInterface = (size - size % 4) / 4;
    }

    // Determine sizes
    label nDomainInterface = nInterfaceTotal - nBaffleInterface;

    // Maximum dimension
    List<int> mapData(max(2 * nInterfaceTotal, 4 * nBaffleInterface), -1);

    bafInterfaces_.setSize(nBaffleInterface);
    domInterfaces_.setSize(nDomainInterface);

    // Face number mapping
    label maxId = max(origFaceId_);
    labelList toFoamFaces(invert(maxId+1, origFaceId_));

    if (nBaffleInterface > 0)
    {
        mapData = -1;

        CCMIOReadOpt1i
        (
            &(globalState_->error),
            interfaceNode,
            "ProstarBaffles",
            mapData.data(),
            kCCMIOStart,
            kCCMIOEnd
        );
        assertNoError("problem reading interface 'ProstarBaffles'");


        // Copy/compress the desired entries (cannot use a SubList)
        // Transform
        //   from [ prostarId faceId1 faceId2 celltableId ]
        //     to [ faceId1 faceId2 ]
        for (label i=0; i < nBaffleInterface; ++i)
        {
            mapData[i*2]   = mapData[i*4+1];
            mapData[i*2+1] = mapData[i*4+2];
        }

        for (label i=2*nBaffleInterface; i < 4*nBaffleInterface; ++i)
        {
            mapData[i] = -1;
        }

        // Translate ccm faceId -> foam faceId
        inplaceRenumber(toFoamFaces, mapData);

        forAll(bafInterfaces_, i)
        {
            bafInterfaces_[i][0] = mapData[i*2];
            bafInterfaces_[i][1] = mapData[i*2+1];
        }
    }

    // Baffles are not domInterfaces, use hash to skip them
    SubList<label> subMap(mapData, 2*nBaffleInterface);
    labelHashSet hashedFace(subMap);

    mapData = -1;

    CCMIOReadOpt2i
    (
        &(globalState_->error),
        interfaceNode,
        "FaceIds",
        mapData.data(),
        kCCMIOStart,
        kCCMIOEnd
    );
    assertNoError("problem reading interface 'FaceIds'");

    // Translate ccm faceId -> foam faceId
    // newer version handles negatives indices
    inplaceRenumber(toFoamFaces, mapData);

    nDomainInterface = 0;
    for (label i=0; i < nInterfaceTotal; ++i)
    {
        label face0 = mapData[i*2];
        label face1 = mapData[i*2+1];

        if
        (
            !hashedFace.found(face0)
         && !hashedFace.found(face1)
        )
        {
            domInterfaces_[nDomainInterface][0] = face0;
            domInterfaces_[nDomainInterface][1] = face1;
            ++nDomainInterface;
            if (nDomainInterface >= domInterfaces_.size())
            {
                break;
            }
        }
    }

    // Truncate for extra safety
    domInterfaces_.setSize(nDomainInterface);
}


// Read monitoring faces
//
void Foam::ccm::reader::readMonitoring
(
    const ccmID& topoId
)
{
#ifdef DEBUG_CCMIOREAD
    Info<< "readMonitoring()" << endl;
#endif

#ifdef WITH_MONITORING
    CCMIONode topoNode, monitorParent;
    CCMIONode monitorNode;

    // CCMIOID -> CCMIONODE
    if
    (
        CCMIOGetEntityNode
        (
            nullptr,
            topoId,
            &topoNode
        )
     == kCCMIONoErr

     // Get "/Meshes/FaceBasedTopology/MonitorBoundaryRegions"
     && CCMIOGetNode
        (
            nullptr,
            topoNode,
            "MonitorBoundaryRegions",
            &monitorParent
        )
     == kCCMIONoErr
    )
    {
        labelList toFoamFaces(invert(nFaces_+1, origFaceId_));

        // Simulate CCMIONextEntity
        for
        (
            int index = 0;
            CCMIOGetNextChildWithLabel
            (
                nullptr,
                monitorParent,
                "boundaryFaces",
                &index,
                &monitorNode
            ) == kCCMIONoErr;
         /* nop */
        )
        {
#ifdef DEBUG_MONITORING
            Info<< "index = " << index << endl;
#endif

            int nMonFaces = 0;
            // Simulate EntitySize ?
            CCMIOReadNodei
            (
                nullptr,
                monitorNode,
                "NumFaces",
                &nMonFaces
            );

            int ccmRegionId = ccmGetEntityIndex(monitorNode);

            ccmID mapId;

            int idVal;
            CCMIOReadNodei
            (
                nullptr,
                monitorNode,
                "MapId",
                &idVal
            );

#ifdef DEBUG_MONITORING
            Info<< "monitoring mapId " << idVal
                << " with nFaces = " << nMonFaces
                << endl;
#endif
            // Is it risky changing parents here?
            CCMIOGetEntity
            (
                nullptr,
                (globalState_->root),
                kCCMIOMap,
                idVal,
                &mapId
            );

            List<label> mapData(nMonFaces);
            // mapData.setSize(nMonFaces);
            // faceCells.setSize(nMonFaces);

            readMap
            (
                mapId,
                mapData
            );

#ifdef DEBUG_MONITORING
            Info<< "map: " << mapData << nl
                << "toFoam: " << toFoamFaces
                << endl;
#endif

            // Translate ccm faceId -> foam faceId
            // newer version handles negatives indices
            inplaceRenumber(toFoamFaces, mapData);

#ifdef DEBUG_MONITORING
            Info<< "map: " << mapData << nl
                << "ccmRegionId: " << ccmRegionId << endl;
#endif

            auto iter = boundaryRegion_.cfind(ccmRegionId);

            word zoneName;
            if (iter.good())
            {
                iter().readEntry("Label", zoneName);
            }
            else
            {
                zoneName = "monitoring_" + Foam::name(ccmRegionId);
            }

            monitoringSets_.insert(zoneName, mapData);

            // CCMIONode subNode;
            //
            //- simulate ReadFaceCells with kCCMIOBoundaryFaces
            // CCMIOGetNode(nullptr, childNode, "Cells", &subNode);
            // CCMIORead1i
            // (
            //     nullptr, subNode, faceCells.data(),
            //     kCCMIOStart, kCCMIOEnd
            // );
            //
            // Info << "cells: " << faceCells << endl;
        }
    }
#endif
}


// Move solid faces from Default_Boundary_Region -> Default_Boundary_Solid
void Foam::ccm::reader::juggleSolids()
{
    if (!option().keepSolid())
    {
        return;
    }

    // Find "Default_Boundary_Region"
    label defaultBoundaryRegion = boundaryRegion_.findIndex
    (
        defaultBoundaryName
    );

    // Find "Default_Boundary_Solid"
    label defaultBoundarySolid = boundaryRegion_.findIndex
    (
        defaultSolidBoundaryName
    );

    // Cannot do anything if Default_Boundary_Region does not exist
    // or if Default_Boundary_Solid already exists
    if
    (
        defaultBoundaryRegion < 0
     || defaultBoundarySolid >= 0
    )
    {
        return;
    }

    // Identify solid cells
    // ~~~~~~~~~~~~~~~~~~~~
    label nSolids = 0;
    bitSet solidCells(cellTableId_.size(), false);
    {
        Map<word> solidMap = cellTable_.solids();

        forAll(cellTableId_, cellI)
        {
            if (solidMap.found(cellTableId_[cellI]))
            {
                solidCells.set(cellI);
                ++nSolids;
            }
        }
    }

    if (!nSolids)
    {
        return;
    }


    // The corresponding Foam patch
    const label patchIndex  = origBndId_.find(defaultBoundaryRegion);
    const label nPatchFaces = patchSizes_[patchIndex];

    labelList patchStarts(patchStartList(nInternalFaces_));
    label adjustPatch = 0;
    for (label i = 0; i < nPatchFaces; ++i)
    {
        label faceI = patchStarts[patchIndex] + i;
        label cellI = faceOwner_[faceI];

        if (solidCells.test(cellI))
        {
            ++adjustPatch;
        }
    }

    // No solid cells on the Default_Boundary_Region
    if (!adjustPatch)
    {
        return;
    }


    // Insert Default_Boundary_Solid immediately after Default_Boundary_Region
    // then we only need to adjust a single patch and can easily re-merge
    // later
    label nPatches = patchSizes_.size();
    patchStarts.setSize(nPatches+1, 0);
    patchSizes_.setSize(nPatches+1, 0);
    origBndId_.setSize(nPatches+1, 0);

    // make room for new entry
    for (label i = nPatches; i > patchIndex; --i)
    {
        patchStarts[i]   = patchStarts[i-1];
        patchSizes_[i]   = patchSizes_[i-1];
        origBndId_[i] = origBndId_[i-1];
    }

    // Adjust start and sizes
    patchSizes_[patchIndex] -= adjustPatch;
    patchSizes_[patchIndex+1] = adjustPatch;
    patchStarts[patchIndex+1] =
        patchStarts[patchIndex] + patchSizes_[patchIndex];

    origBndId_[patchIndex+1] = boundaryRegion_.append
    (
        dictionary
        (
            IStringStream
            (
                "BoundaryType wall;"
                "Label " + word(defaultSolidBoundaryName) + ";"
            )()
        )
    );

    label fluidFace = patchStarts[patchIndex];
    label solidFace = patchStarts[patchIndex+1];

    labelList oldToNew(identity(nFaces_));
    for (label i = 0; i < nPatchFaces; ++i)
    {
        label faceI = patchStarts[patchIndex] + i;
        label cellI = faceOwner_[faceI];

        if (solidCells.test(cellI))
        {
            oldToNew[faceI] = solidFace++;
        }
        else
        {
            oldToNew[faceI] = fluidFace++;
        }
    }

    // Re-order faces, owners/neighbours
    inplaceReorder(oldToNew, faces_);
    inplaceReorder(oldToNew, faceOwner_);
    inplaceReorder(oldToNew, faceNeighbour_);
    inplaceReorder(oldToNew, origFaceId_);

    renumberInterfaces(oldToNew);
}


// In CCM, separate mesh domains are used for fluid/porous/solid
// Thus the fluid/porous/solid cells correspond uniquely to a face owner
void Foam::ccm::reader::removeUnwanted()
{
    // Identify fluid/porous/solid cells for removal
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    label nRemove = 0;
    bitSet removeCells(cellTableId_.size(), false);

    {
        Map<word> fluidMap  = cellTable_.fluids();
        Map<word> porousMap = selectPorous(cellTable_);
        Map<word> solidMap  = cellTable_.solids();
        Map<word> removeMap;

        forAll(cellTableId_, cellI)
        {
            label tableId = cellTableId_[cellI];

            if
            (
                porousMap.found(tableId)
              ? !option().keepPorous()
              : fluidMap.found(tableId)
              ? !option().keepFluid()
              : solidMap.found(tableId)
              ? !option().keepSolid()
              : false
            )
            {
                removeCells.set(cellI);
                ++nRemove;
                removeMap.set(tableId, cellTable_.name(tableId));
            }
        }

        if (nRemove)
        {
            Map<word> keepMap;

            forAllConstIters(cellTable_, iter)
            {
                const label tableId = iter.key();
                if (!removeMap.found(tableId))
                {
                    keepMap.set(tableId, cellTable_.name(tableId));
                }
            }

            Info<<"remove "<< nRemove << " cells in "
                << removeMap.size() << " unwanted cellZone(s)" << nl;

            forAllConstIters(removeMap, iter)
            {
                Info<< "    zone "
                    << iter.key() << " : " << iter.val() << nl;
            }

            Info<<"retain "<< (nCells_ - nRemove) << " cells in "
                << keepMap.size() << " cellZone(s)" << nl;

            forAllConstIters(keepMap, iter)
            {
                Info<< "    zone "
                    << iter.key() << " : " << iter.val() << nl;
            }
        }
    }

    if (!nRemove)
    {
        return;
    }

    // Remove all faces where the owner corresponds to a removed cell
    // Adjust the nInternalFaces and patch sizes accordingly
    label adjustInternal = 0;
    labelList adjustPatchSize(patchSizes_.size(), Zero);

    label newFaceI = 0;
    label oldFaceI = nFaces_ - 1;
    labelList oldToNew(nFaces_, -1);
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        label cellI = faceOwner_[faceI];
        if (removeCells.test(cellI))
        {
            if (faceI < nInternalFaces_)
            {
                ++adjustInternal;
            }
            else
            {
                // need to adjust the patch sizes
                label beg = nInternalFaces_;
                forAll(patchSizes_, patchI)
                {
                    label end = beg + patchSizes_[patchI];

                    if (faceI >= beg && faceI < end)
                    {
                        ++adjustPatchSize[patchI];
                        break;
                    }

                    beg = end;
                }
            }

            // Put discarded faces at the end of the list
            oldToNew[faceI] = oldFaceI--;
        }
        else
        {
            if (newFaceI != faceI)
            {
                faces_[newFaceI] = faces_[faceI];
                faceOwner_[newFaceI] = faceOwner_[faceI];
                faceNeighbour_[newFaceI] = faceNeighbour_[faceI];
                origFaceId_[newFaceI] = origFaceId_[faceI];
            }

            // New position for the face
            oldToNew[faceI] = newFaceI++;
        }
    }

    nFaces_ = newFaceI;

    // Redimension
    faces_.setSize(nFaces_);
    faceOwner_.setSize(nFaces_);
    faceNeighbour_.setSize(nFaces_);
    origFaceId_.setSize(nFaces_);

    // Adjust internal faces and patch sizes

    nInternalFaces_ -= adjustInternal;
    forAll(patchSizes_, patchI)
    {
        patchSizes_[patchI] -= adjustPatchSize[patchI];
    }

    renumberInterfaces(oldToNew);

    // Renumber cells
    label nCell = 0;
    oldToNew.setSize(nCells_, -1);
    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        if (!removeCells.test(cellI))
        {
            if (nCell != cellI)
            {
                origCellId_[nCell] = origCellId_[cellI];
                cellTableId_[nCell] = cellTableId_[cellI];
            }
            oldToNew[cellI] = nCell;
            ++nCell;
        }
    }

    inplaceRenumber(oldToNew, faceOwner_);
    inplaceRenumber(oldToNew, faceNeighbour_);

    // Redimension
    nCells_ = nCell;
    origCellId_.setSize(nCells_);
    cellTableId_.setSize(nCells_);


    // Remove unused points - adjust points, faces accordingly
    oldToNew.setSize(nPoints_);
    oldToNew = -1;

    // Mark up all the used points
    for (const face& f : faces_)
    {
        for (const label pointi : f)
        {
            ++oldToNew[pointi];
        }
    }

    label nPointUsed = 0;
    forAll(oldToNew, ptI)
    {
        if (oldToNew[ptI] >= 0)
        {
            oldToNew[ptI] = nPointUsed;
            if (ptI != nPointUsed)
            {
                points_[nPointUsed] = points_[ptI];
            }
            ++nPointUsed;
        }
    }

    nPoints_ = nPointUsed;
    points_.setSize(nPoints_);

    for (face& f : faces_)
    {
        inplaceRenumber(oldToNew, f);
    }

    // Report sizes
    printSizes();
}


void Foam::ccm::reader::validateInterface
(
    List<labelPair>& lst
)
{
    label nElem = 0;
    forAll(lst, elemI)
    {
        label face0 = lst[elemI][0];
        label face1 = lst[elemI][1];

        if (face0 < nFaces_ && face1 < nFaces_)
        {
            if (nElem != elemI)
            {
                lst[nElem][0] = face0;
                lst[nElem][1] = face1;
            }
            ++nElem;
        }
    }
    lst.setSize(nElem);
}


void Foam::ccm::reader::renumberInterfaces
(
    const labelUList& oldToNew
)
{
    forAll(domInterfaces_, elemI)
    {
        domInterfaces_[elemI][0] = oldToNew[domInterfaces_[elemI][0]];
        domInterfaces_[elemI][1] = oldToNew[domInterfaces_[elemI][1]];
    }

    forAll(bafInterfaces_, elemI)
    {
        bafInterfaces_[elemI][0] = oldToNew[bafInterfaces_[elemI][0]];
        bafInterfaces_[elemI][1] = oldToNew[bafInterfaces_[elemI][1]];
    }

    validateInterface(domInterfaces_);
    validateInterface(bafInterfaces_);
}


//
// 1) remove interfaces between domains (fluid/porosity; fluid/solid, etc)
// 2) reorganize baffle interfaces into [0-N/2; N/2-N] lists at the beginning
//    of the corresponding patch
//
void Foam::ccm::reader::cleanupInterfaces()
{
    validateInterface(bafInterfaces_);
    validateInterface(domInterfaces_);

    if (bafInterfaces_.size() <= 0 && domInterfaces_.size() <= 0)
    {
        Info<<"0 baffle interface pairs" << nl
            <<"0 domain interface pairs" << endl;
        return;
    }

#ifdef DEBUG_BAFFLES
    Info<< "baffle Interfaces " << bafInterfaces_ << nl
        << "domain Interfaces " << domInterfaces_ << nl
        << "nCells:" << nCells_ << nl
        << "nFaces:" << nFaces_ << nl
        << "patchSizes:"  << patchSizes_ << nl
        << "nInternalFaces:" << nInternalFaces_ << endl;

    forAll(domInterfaces_, elemI)
    {
        const label face0 = domInterfaces_[elemI][0];
        const label face1 = domInterfaces_[elemI][1];

        Info<< "interface [" << elemI << "] = "
            << face0 << " - " << face1 << " own/neigh = "
            << faceOwner_[face0] << "/" << faceNeighbour_[face0] << "  "
            << faceOwner_[face1] << "/" << faceNeighbour_[face1] << endl;
    }
#endif

    // Only reorder faces that need it
    labelList oldToNew(nFaces_, -1);

    // - move side0 face from domInterfaces to join the internal faces
    // - move the redundant side1 to the end of the list for later deletion
    label begOfList = nInternalFaces_;
    label endOfList = nFaces_ - 1;

    // The patch sizes (and the start) will definitely change
    const labelList origPatchStarts(patchStartList(nInternalFaces_));
    labelList adjustPatchSize(patchSizes_.size(), Zero);
    labelList bafflePatchCount(patchSizes_.size(), Zero);

    // The new dimensions after merging the domain interfaces:
    nInternalFaces_ += domInterfaces_.size();
    nFaces_         -= domInterfaces_.size();

    Info<< domInterfaces_.size() << " domain interface pairs";
    if (domInterfaces_.size())
    {
        Info<<" to merge" << endl;
        printSizes();
    }
    else
    {
        Info<< endl;
    }

    forAll(domInterfaces_, elemI)
    {
        label face0 = domInterfaces_[elemI][0];
        label face1 = domInterfaces_[elemI][1];

        oldToNew[face0] = begOfList++;
        oldToNew[face1] = endOfList--;     // End of list for truncationx

        // face0 gets a new neighbour, face1 loses its owner
        faceNeighbour_[face0] = faceOwner_[face1];
        faceOwner_[face1] = -1;

        // Need to adjust the patch sizes
        forAll(patchSizes_, patchI)
        {
            label beg = origPatchStarts[patchI];
            label end = beg + patchSizes_[patchI];

            if (face0 >= beg && face0 < end)
            {
                ++adjustPatchSize[patchI];
            }
            if (face1 >= beg && face1 < end)
            {
                ++adjustPatchSize[patchI];
            }
        }
    }

    // Count the number of baffles per patch
    forAll(bafInterfaces_, elemI)
    {
        label face0 = bafInterfaces_[elemI][0];
        label face1 = bafInterfaces_[elemI][1];

        forAll(patchSizes_, patchI)
        {
            label beg = origPatchStarts[patchI];
            label end = beg + patchSizes_[patchI];

            if (face0 >= beg && face0 < end)
            {
                ++bafflePatchCount[patchI];
            }
            if (face1 >= beg && face1 < end)
            {
                ++bafflePatchCount[patchI];
            }
        }
    }


    if (option().removeBaffles())
    {
        // The new dimensions after merging the baffles:
        nInternalFaces_ += bafInterfaces_.size();
        nFaces_         -= bafInterfaces_.size();
    }

    Info<< bafInterfaces_.size() << " baffle interface pairs";
    if (bafInterfaces_.size())
    {
        if (option().removeBaffles())
        {
            Info<< " to merge" << endl;
            printSizes();
        }
        else
        {
            Info<< " to be sorted" << endl;
        }
    }
    else
    {
        Info<< endl;
    }

    if (option().removeBaffles())
    {
        forAll(bafInterfaces_, elemI)
        {
            label face0 = bafInterfaces_[elemI][0];
            label face1 = bafInterfaces_[elemI][1];

            oldToNew[face0] = begOfList++;
            oldToNew[face1] = endOfList--;     // End of list for truncation

            // face0 gets a new neighbour, face1 loses its owner
            faceNeighbour_[face0] = faceOwner_[face1];
            faceOwner_[face1] = -1;
        }


        // Boundary faces continue from the new nInternalFaces
        label pos = nInternalFaces_;
        forAll(patchSizes_, patchI)
        {
            label beg = origPatchStarts[patchI];
            label end = beg + patchSizes_[patchI];

            // Fill in values for the remainder of the boundary patch
            for (label faceI = beg; faceI < end; ++faceI)
            {
                if (oldToNew[faceI] < 0)
                {
                    oldToNew[faceI] = pos;
                    ++pos;
                }
            }
        }

        // Baffles have been resolved - remove last traces
        bafInterfaces_.clear();
    }
    else
    {
        // This check is probably unnecessary
        forAll(bafflePatchCount, patchI)
        {
            if (bafflePatchCount[patchI] % 2)
            {
                Info<< "WARNING: patch " << patchI
                    << " has an uneven number of baffles ("
                    << bafflePatchCount[patchI] << ") expect strange results"
                    << endl;
            }
        }


        // Reordered faces continue from the new nInternalFaces
        label pos = nInternalFaces_;
        forAll(patchSizes_, patchI)
        {
            const label beg = origPatchStarts[patchI];
            const label end = beg + patchSizes_[patchI];

            const label nsize = bafflePatchCount[patchI];
            if (nsize > 0)
            {
                // Reorganize baffle interfaces into [0-N/2; N/2-N] lists
                // at the beginning of the corresponding patch
                const label nsizeby2 = (nsize - nsize % 2) / 2;
                label nsorted = 0;

                // Renumber the normal (baffle) interfaces
                forAll(bafInterfaces_, elemI)
                {
                    const label face0 = bafInterfaces_[elemI][0];
                    const label face1 = bafInterfaces_[elemI][1];

                    if
                    (
                        (face0 >= beg && face0 < end)
                     || (face1 >= beg && face1 < end)
                    )
                    {
                        oldToNew[face0] = pos + nsorted;
                        oldToNew[face1] = pos + nsorted + nsizeby2;

                        // Mark destination of the faces, but cannot renumber
                        // yet. Use negative to potential overlap with other
                        // patch regions
                        bafInterfaces_[elemI][0] = -oldToNew[face0];
                        bafInterfaces_[elemI][1] = -oldToNew[face1];

                        ++nsorted;
                    }
                }
                pos += 2*nsorted;
            }

            // Fill in values for the remainder of the boundary patch
            for (label faceI = beg; faceI < end; ++faceI)
            {
                if (oldToNew[faceI] < 0)
                {
                    oldToNew[faceI] = pos;
                    ++pos;
                }
            }
        }

        // Finalize new numbers for the normal (baffle) interfaces
        forAll(bafInterfaces_, elemI)
        {
            bafInterfaces_[elemI][0] = abs(bafInterfaces_[elemI][0]);
            bafInterfaces_[elemI][1] = abs(bafInterfaces_[elemI][1]);
        }
    }

#ifdef DEBUG_BAFFLES
    Info<< "remap with " << oldToNew << nl
        << "owners:" << faceOwner_ << nl
        << "neighbours:" << faceNeighbour_ << nl
        << endl;
#endif

    // Re-order faces, owners/neighbours
    inplaceReorder(oldToNew, faces_);
    inplaceReorder(oldToNew, faceOwner_);
    inplaceReorder(oldToNew, faceNeighbour_);
    inplaceReorder(oldToNew, origFaceId_);

    if (monitoringSets_.size())
    {
#ifdef WITH_MONITORING
        // Modify oldToNew mapping to account for monitoring faces that
        // coincided with a domain interface
        //
        // TODO - should modify flip map as well
        forAll(domInterfaces_, elemI)
        {
            label face0 = domInterfaces_[elemI][0];
            label face1 = domInterfaces_[elemI][1];
            oldToNew[face1] = oldToNew[face0];
        }

        forAllIters(monitoringSets_, iter)
        {
            inplaceRenumber(oldToNew, iter.val());
        }
#endif
    }

    // We can finally drop this information now
    domInterfaces_.clear();

    // Truncate lists
    faces_.setSize(nFaces_);
    faceOwner_.setSize(nFaces_);
    faceNeighbour_.setSize(nFaces_);
    origFaceId_.setSize(nFaces_);

    // Remove empty patches:

    // Fix patch sizes:
    oldToNew.setSize(patchSizes_.size());
    oldToNew = -1;

    label nPatches = 0;
    forAll(patchSizes_, patchI)
    {
        patchSizes_[patchI] -= adjustPatchSize[patchI];
        if (option().removeBaffles())
        {
            patchSizes_[patchI] -= bafflePatchCount[patchI];
        }

        if (patchSizes_[patchI])
        {
            oldToNew[patchI] = nPatches++;
        }
    }

    inplaceReorder(oldToNew, patchSizes_);
    inplaceReorder(oldToNew, origBndId_);

    patchSizes_.setSize(nPatches);
    origBndId_.setSize(nPatches);

#ifdef DEBUG_BAFFLES
    Info<< "nCells:" << nCells_ << nl
        << "nFaces:" << nFaces_ << nl
        << "PatchSizes:"  << patchSizes_ << nl
        << "nInternalFaces:" << nInternalFaces_ << nl
        << endl;
#endif
}


//
// Merge STARCCM in-place interfaces
//
void Foam::ccm::reader::mergeInplaceInterfaces()
{
    if (interfaceDefinitions_.empty())
    {
        return;
    }
    if (!option().mergeInterfaces())
    {
        Info<< interfaceDefinitions_.size() << " interface definitions"
            << " - leaving unmerged" << endl;
        return;
    }

    // List of patch pairs that are interfaces
    DynamicList<labelPair> interfacePatches(interfaceDefinitions_.size());

    label nWarn = 0;

    forAllConstIters(interfaceDefinitions_, iter)
    {
        const interfaceEntry& ifentry = iter.val();

        labelPair patchPair
        (
            origBndId_.find(ifentry.bnd0),
            origBndId_.find(ifentry.bnd1)
        );

        if
        (
            patchPair[0] == patchPair[1]
         || patchPair[0] < 0
         || patchPair[1] < 0
        )
        {
            // This should not happen
            Info<<"Warning : bad interface " << ifentry.id << " " << ifentry
                <<" on patches " << patchPair << endl;
        }
        else if
        (
            patchSizes_[patchPair[0]] != patchSizes_[patchPair[1]]
         || patchSizes_[patchPair[0]] == 0
         || patchSizes_[patchPair[1]] == 0
        )
        {
            if (!nWarn++)
            {
                Info<<"Warning: skip interface with zero or different"
                    << " number of faces" << nl;
            }

            Info<<"  Interface:" << ifentry.id << " " << ifentry
                <<" patches " << patchPair
                <<" sizes ("
                << patchSizes_[patchPair[0]]
                << " " << patchSizes_[patchPair[1]] << ")"
                << nl;
        }
        else
        {
            interfacePatches.append(patchPair);
        }
    }

    if (interfacePatches.empty())
    {
        return;
    }


    // Local point mapping
    labelList mergedPointMap;

    // Global remapping
    labelList oldToNew(identity(points_.size()));

    const labelList origPatchStarts(patchStartList(nInternalFaces_));

    label nMergedTotal = 0;

    // Markup points to merge
    bitSet whichPoints(points_.size());

    Info<< "interface merge points (tol="
        << option().mergeTol() << "):" << endl;

    DynamicList<label> interfacesToMerge(interfacePatches.size());
    forAll(interfacePatches, interI)
    {
        const label patch0 = interfacePatches[interI][0];
        const label patch1 = interfacePatches[interI][1];
        const label nPatch0Faces = patchSizes_[patch0];
        const label nPatch1Faces = patchSizes_[patch1];

        // Markup points to merge
        whichPoints.reset();
        for (label local0FaceI = 0; local0FaceI < nPatch0Faces; ++local0FaceI)
        {
            const face& f = faces_[origPatchStarts[patch0] + local0FaceI];

            for (const label pointi : f)
            {
                // Simultaneously account for previous point merges
                whichPoints.set(oldToNew[pointi]);
            }
        }
        for (label local1FaceI = 0; local1FaceI < nPatch1Faces; ++local1FaceI)
        {
            const face& f = faces_[origPatchStarts[patch1] + local1FaceI];

            for (const label pointi : f)
            {
                // Simultaneously account for previous point merges
                whichPoints.set(oldToNew[pointi]);
            }
        }

        // The global addresses
        labelList addr(whichPoints.toc());

        const UIndirectList<point> pointsToMerge(points_, addr);

        Info<< "    patch "  << patch0 << ',' << patch1 << ": ("
            << nPatch0Faces << " and " << nPatch1Faces << " faces) " << flush;

        const label nMerged = mergePoints
        (
            pointsToMerge,
            option().mergeTol(),
            false,
            mergedPointMap
        );

        Info<< nMerged << " from " << pointsToMerge.size() << " points"
            << endl;

        if (nMerged)
        {
            // Two-steps:
            // * Identify duplicate points - without changing the order!
            // * Transcribe local to global addressing (oldToNew)

            forAll(mergedPointMap, pti)
            {
                const label mergedPti = mergedPointMap[pti];

                if (mergedPti < 0)
                {
                    continue;  // Already seen as duplicate
                }

                const label origPointi = oldToNew[addr[pti]];

                // Find further duplicate points
                for
                (
                    label dupPti = pti+1;
                    (dupPti = mergedPointMap.find(mergedPti, dupPti)) != -1;
                    ++dupPti
                )
                {
                    oldToNew[addr[dupPti]] = origPointi;
                    mergedPointMap[dupPti] = -1;  // Mark as already seen
                }
            }

            interfacesToMerge.append(interI);
            nMergedTotal += nMerged;
        }
    }


    //
    // Nothing to do
    //
    if (!nMergedTotal)
    {
        return;
    }

    // Update point references to account for point merge:
    for (face& f : faces_)
    {
        inplaceRenumber(oldToNew, f);
    }

    // Determine which points are actually in use:
    oldToNew.resize_nocopy(nPoints_);
    oldToNew = -1;

    // Mark up all the used points
    for (const face& f : faces_)
    {
        for (const label pointi : f)
        {
            ++oldToNew[pointi];
        }
    }

    label nPointUsed = 0;
    forAll(oldToNew, ptI)
    {
        if (oldToNew[ptI] >= 0)
        {
            oldToNew[ptI] = nPointUsed;
            if (ptI != nPointUsed)
            {
                points_[nPointUsed] = points_[ptI];
            }
            ++nPointUsed;
        }
    }

    // Info<< "merge " << nMergedTotal << " points from "
    //     << nPoints_ << " to " << nPointUsed << endl;

    nPoints_ = nPointUsed;
    points_.resize(nPoints_);

    for (face& f : faces_)
    {
        inplaceRenumber(oldToNew, f);
    }


    //
    // Merge the faces as well
    //
    Info<< "interface merge faces:" << endl;

    nMergedTotal = 0;
    labelList adjustPatchSize(patchSizes_.size(), Zero);
    forAll(interfacesToMerge, mergeI)
    {
        const label patch0 = interfacePatches[interfacesToMerge[mergeI]][0];
        const label patch1 = interfacePatches[interfacesToMerge[mergeI]][1];

        labelList faceAddr0(patchSizes_[patch0]);
        labelList faceAddr1(patchSizes_[patch1]);

        forAll(faceAddr0, localFaceI)
        {
            faceAddr0[localFaceI] = origPatchStarts[patch0] + localFaceI;
        }
        forAll(faceAddr1, localFaceI)
        {
            faceAddr1[localFaceI] = origPatchStarts[patch1] + localFaceI;
        }

        if (faceAddr0.size() != faceAddr1.size())
        {
            // This should not occur, we avoided the same thing above
            continue;
        }

        // Improve comparison speed by sorting by distance
        SortableList<scalar> pts0MagSqr
        (
            magSqr
            (
                uindirectPrimitivePatch
                (
                    UIndirectList<face>
                    (
                        faces_,
                        faceAddr0
                    ),
                    points_
                ).faceCentres()
            )
        );
        SortableList<scalar> pts1MagSqr
        (
            magSqr
            (
                uindirectPrimitivePatch
                (
                    UIndirectList<face>
                    (
                        faces_,
                        faceAddr1
                    ),
                    points_
                ).faceCentres()
            )
        );

        label nMerged = 0;

        // Record which faces failed to merge - use slower ad hoc merging
        labelHashSet failed0, failed1;
        forAll(pts0MagSqr, sortI)
        {
            const label face0I = faceAddr0[pts0MagSqr.indices()[sortI]];
            const label face1I = faceAddr1[pts1MagSqr.indices()[sortI]];

            // This is what we expect
            if (face::compare(faces_[face0I], faces_[face1I]))
            {
                ++nMerged;

                // Moved from boundary patch to internal patch
                ++adjustPatchSize[patch0];
                ++adjustPatchSize[patch1];

                if (faceOwner_[face0I] < faceOwner_[face1I])
                {
                    // keep 0, discard 1
                    faceNeighbour_[face0I] = faceOwner_[face1I];
                    faceNeighbour_[face1I] = faceOwner_[face1I] = -1;
                }
                else
                {
                    // keep 1, discard 0
                    faceNeighbour_[face1I] = faceOwner_[face0I];
                    faceNeighbour_[face0I] = faceOwner_[face0I] = -1;
                }
            }
            else
            {
                failed0.set(face0I);
                failed1.set(face1I);
            }
        }

        // Perhaps some sorting issues, recheck
        if (failed0.size())
        {
            // Note which one were successful
            labelHashSet done(failed0.size());

            for (const label face0I : failed0)
            {
                for (const label face1I : failed1)
                {
                    // This is what we expect
                    if (face::compare(faces_[face0I], faces_[face1I]))
                    {
                        ++nMerged;

                        // Moved from boundary patch to internal patch
                        ++adjustPatchSize[patch0];
                        ++adjustPatchSize[patch1];

                        if (faceOwner_[face0I] < faceOwner_[face1I])
                        {
                            // keep 0, discard 1
                            faceNeighbour_[face0I] = faceOwner_[face1I];
                            faceNeighbour_[face1I] = faceOwner_[face1I] = -1;
                        }
                        else
                        {
                            // keep 1, discard 0
                            faceNeighbour_[face1I] = faceOwner_[face0I];
                            faceNeighbour_[face0I] = faceOwner_[face0I] = -1;
                        }

                        failed1.erase(face1I);   // Never check again
                        done.set(face0I);        // Mark as done
                        break;                   // Stop looking
                    }
                }
            }

            // Transfer to note how many were successful
            failed0 = done;
        }

        Info<< "    patch "  << patch0 << ',' << patch1 << ": "
            << nMerged << " from " << faceAddr0.size() << " faces";

        if (failed0.size())
        {
            Info<< " (" << failed0.size() << " merged ad hoc)";
        }
        Info<< endl;


        nMergedTotal += nMerged;
    }


    // Nothing to do
    if (!nMergedTotal)
    {
        return;
    }

    // Info<< "merge " << nMergedTotal << " faces from "
    //     << nFaces_ << " to " << (nFaces_ - nMergedTotal) << endl;

    oldToNew.setSize(nFaces_);
    oldToNew = -1;

    // Remaining external faces will be shifted here:
    label extFaceI = nInternalFaces_ + nMergedTotal;

    // Start over
    nInternalFaces_ = 0;
    label nFaceUsed = 0;
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        if (faceOwner_[faceI] != -1)
        {
            if (faceNeighbour_[faceI] != -1)
            {
                // Internal face
                oldToNew[faceI] = nInternalFaces_;
                ++nInternalFaces_;
                ++nFaceUsed;
            }
            else
            {
                // External face
                oldToNew[faceI] = extFaceI;
                ++extFaceI;
                ++nFaceUsed;
            }
        }
    }

    if (nFaceUsed != extFaceI)
    {
        FatalErrorInFunction
            << "coding error: used " << nFaceUsed
            << " faces, but expected to use " << extFaceI << " faces"
            << exit(FatalError);
    }

    // Re-order faces, owners/neighbours
    inplaceReorder(oldToNew, faces_, true);
    inplaceReorder(oldToNew, faceOwner_, true);
    inplaceReorder(oldToNew, faceNeighbour_, true);
    inplaceReorder(oldToNew, origFaceId_, true);

    nFaces_ = nFaceUsed;

    faces_.setSize(nFaces_);
    faceOwner_.setSize(nFaces_);
    faceNeighbour_.setSize(nFaces_);
    origFaceId_.setSize(nFaces_);

    // Fix patch sizes:
    oldToNew.setSize(patchSizes_.size());
    oldToNew = -1;

    label nPatches = 0;
    forAll(patchSizes_, patchI)
    {
        patchSizes_[patchI] -= adjustPatchSize[patchI];
        if (patchSizes_[patchI])
        {
            oldToNew[patchI] = nPatches++;
        }
    }

    inplaceReorder(oldToNew, patchSizes_);
    inplaceReorder(oldToNew, origBndId_);

    patchSizes_.setSize(nPatches);
    origBndId_.setSize(nPatches);

    // Report we are done
    Info<< ".." << endl;
}


//
// Re-order mesh into Foam convention
// - owner < neighbour
// - face vertices such that normal points away from owner
// - order faces: upper-triangular for internal faces;
//    boundary faces after internal faces
//
void Foam::ccm::reader::reorderMesh()
{
    // Set owner/neighbour so owner < neighbour
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(faceOwner_, faceI)
    {
        const label nbr = faceNeighbour_[faceI];
        const label own = faceOwner_[faceI];

        if (nbr >= cellTableId_.size() || own >= cellTableId_.size())
        {
            FatalErrorInFunction
                << "face:" << faceI
                << " nbr:" << nbr
                << " own:" << own
                << " nCells:" << cellTableId_.size()
                << exit(FatalError);
        }

        if (nbr >= 0 && nbr < own)
        {
            faceOwner_[faceI] = faceNeighbour_[faceI];
            faceNeighbour_[faceI] = own;
            faces_[faceI].flip();
        }

        // And check the face
        const face& f = faces_[faceI];

        for (const label pointi : f)
        {
            if (pointi < 0 || pointi >= points_.size())
            {
                FatalErrorInFunction
                    << "face:" << faceI << " f:" << f
                    << abort(FatalError);
            }
        }
    }

    // Do upper-triangular ordering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList oldToNew(faceOwner_.size(), -1);

    // Create cells (inverse of face-to-cell addressing)
    cellList cellFaceAddr;
    primitiveMesh::calcCells
    (
        cellFaceAddr,
        faceOwner_,
        faceNeighbour_,
        nCells_
    );

    label newFaceI = 0;
    forAll(cellFaceAddr, cellI)
    {
        const labelList& cFaces = cellFaceAddr[cellI];
        SortableList<label> nbr(cFaces.size(), -1);

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];
            label nbrCellI = faceNeighbour_[faceI];

            if (nbrCellI >= 0)
            {
                // Internal face. Get cell on other side
                if (nbrCellI == cellI)
                {
                    nbrCellI = faceOwner_[faceI];
                }

                // cellI is master
                if (cellI < nbrCellI)
                {
                    nbr[i] = nbrCellI;
                }
            }
        }

        nbr.sort();

        forAll(nbr, i)
        {
            if (nbr[i] >= 0)
            {
                oldToNew[cFaces[nbr.indices()[i]]] = newFaceI++;
            }
        }
    }

    // Reorder faces accordingly
    inplaceReorder(oldToNew, faces_);
    inplaceReorder(oldToNew, faceOwner_);
    inplaceReorder(oldToNew, faceNeighbour_);
    inplaceReorder(oldToNew, origFaceId_);

    forAllIters(monitoringSets_, iter)
    {
        labelList& lst = iter.val();
        inplaceRenumber(oldToNew, lst);

        // disallow monitoring on boundaries
        label nElem = 0;
        forAll(lst, i)
        {
            if (lst[i] >= 0 && lst[i] < nInternalFaces_)
            {
                if (nElem != i)
                {
                    lst[nElem] = lst[i];
                }
                ++nElem;
            }
        }

        if (nElem)
        {
            lst.setSize(nElem);
        }
        else
        {
            Info << "remove monitor " << iter.key() << endl;
            monitoringSets_.erase(iter);
        }
    }
}


// Attach patches
//
// - patchSizes_   : obvious
// - origBndId_ : lookup for patch name/type
//
void Foam::ccm::reader::addPatches
(
    polyMesh& mesh
) const
{
    // Create patches
    // use patch types to determine what Foam types to generate
    polyPatchList newPatches(origBndId_.size());

    label meshFaceI = nInternalFaces_;
    wordHashSet hashedNames(origBndId_.size());

    // lookup patch names/types from the problem description
    // provide some fallback values
    forAll(newPatches, patchI)
    {
        const word fallbackName(polyPatch::defaultName(patchI));
        word patchName;
        word patchType;

        auto citer = boundaryRegion_.cfind(origBndId_[patchI]);

        if (citer.good())
        {
            citer().readEntry("Label", patchName);
            citer().readEntry("BoundaryType", patchType);
        }
        else
        {
            patchName = fallbackName;
            patchType = "patch";
        }

        // Avoid duplicate names
        // - don't bother checking if the modified name is also a duplicate
        if (hashedNames.found(patchName))
        {
            Info<< "renamed patch " << patchName << " to ";
            patchName = fallbackName + "_" + patchName;
            Info<< patchName << endl;
        }
        hashedNames.insert(patchName);

        Info<< "patch " << patchI
            << " (start: " << meshFaceI << " size: " << patchSizes_[patchI]
            << ") name: " << patchName
            << endl;

        if (patchType == "wall")
        {
            newPatches.set
            (
                patchI,
                new wallPolyPatch
                (
                    patchName,
                    patchSizes_[patchI],
                    meshFaceI,
                    patchI,
                    mesh.boundaryMesh(),
                    patchType
                )
            );
        }
        else if (patchType == "symmetry")
        {
            newPatches.set
            (
                patchI,
                new symmetryPolyPatch
                (
                    patchName,
                    patchSizes_[patchI],
                    meshFaceI,
                    patchI,
                    mesh.boundaryMesh(),
                    patchType
                )
            );
        }
        else if (patchType == "empty")
        {
            // Note: not ccm name, may have been introduced by us
            newPatches.set
            (
                patchI,
                new emptyPolyPatch
                (
                    patchName,
                    patchSizes_[patchI],
                    meshFaceI,
                    patchI,
                    mesh.boundaryMesh(),
                    patchType
                )
            );
        }
        else
        {
            // All other ccm types become straight polyPatch:
            // 'inlet', 'outlet', 'pressure'.
            newPatches.set
            (
                patchI,
                new polyPatch
                (
                    patchName,
                    patchSizes_[patchI],
                    meshFaceI,
                    patchI,
                    mesh.boundaryMesh(),
                    patchType
                )
            );
        }

        meshFaceI += patchSizes_[patchI];
    }

    if (meshFaceI != mesh.nFaces())
    {
        FatalErrorInFunction
            << "meshFaceI:" << meshFaceI << " nFaces:" << mesh.nFaces()
            << abort(FatalError);
    }

    mesh.addPatches(newPatches);
}



// Attach faceZones based on the monitoring boundary conditions
void Foam::ccm::reader::addFaceZones
(
    polyMesh& mesh
) const
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

        ++nZone;
    }

    mesh.faceZones().writeOpt(IOobject::AUTO_WRITE);
    warnDuplicates("faceZones", mesh.faceZones().names());
}


// Remove most of the ccm-specific information with the exception of information
// auxiliary to the normal polyMesh:
//   (cellTableId_, origCellId_, origFaceId_, interfaces_, ...)
void Foam::ccm::reader::clearGeom()
{
    // allow re-reading the file
    if (geometryStatus_ == OKAY || geometryStatus_ == READ)
    {
        geometryStatus_ = UNKNOWN;
    }

    points_.clear();
    faces_.clear();
    faceOwner_.clear();
    faceNeighbour_.clear();

    origBndId_.clear();
    patchSizes_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyMesh> Foam::ccm::reader::mesh
(
    const objectRegistry& registry,
    const fileName& remappingDictName
)
{
    if (!readGeometry())
    {
        return nullptr;
    }

    // merge cellTable and rename boundaryRegion
    remapMeshInfo(registry, remappingDictName);

    // Construct polyMesh
    // ~~~~~~~~~~~~~~~~~~
    auto meshPtr = autoPtr<polyMesh>::New
    (
        IOobject
        (
            polyMesh::defaultRegion,
            "constant",
            registry,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        std::move(points_),
        std::move(faces_),
        std::move(faceOwner_),
        std::move(faceNeighbour_)
    );
    polyMesh& mesh = *meshPtr;

    addPatches(mesh);

    // Attach cellZones based on the cellTable Id
    // any other values can be extracted later from the cellTable dictionary
    cellTable_.addCellZones(mesh, cellTableId_);
    warnDuplicates("cellZones", mesh.cellZones().names());

    addFaceZones(mesh);

    warnDuplicates("boundaries", mesh.boundaryMesh().names());
    clearGeom();

    return meshPtr;
}


// ************************************************************************* //
