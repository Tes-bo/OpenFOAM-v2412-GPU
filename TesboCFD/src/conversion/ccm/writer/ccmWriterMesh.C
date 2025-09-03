/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "ccmWriter.H"
#include "ccmInternal.H"  // include last to avoid any strange interactions


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::ccm::writer::prostarCellFaceId
(
    const label& cellId,
    const label& faceI
) const
{
    const faceList&   faces = mesh_.faces();
    const cellShape&  shape = mesh_.cellShapes()[cellId];
    const labelList& cFaces = mesh_.cells()[cellId];

    label cellFaceId = cFaces.find(faceI);
    label mapIndex = shape.model().index();

    if (ccm::debug > 1)
    {
        Info<< " face[" << faceI << "]: " << faces[faceI] << nl
            << " owner: " << cellId
            << " cellFace: " << cellFaceId
            << " cellFaces: " << cFaces
            << " shape [index " << mapIndex << "] "  << shape
            << endl;

        Info << "cellFaces" << nl;
        forAll(cFaces, cFaceI)
        {
            Info<< "  face [" << cFaces[cFaceI] << "] = "
                << faces[cFaces[cFaceI]] << nl;
        }

        Info << "shapeFaces" << nl;
        {
            const faceList sFaces = shape.faces();
            forAll(sFaces, sFaceI)
            {
                Info<< "  sFace[" << sFaceI << "] = "
                    << sFaces[sFaceI] << nl;
            }
        }
    }

    // The face order returned by primitiveMesh::cells() is not
    // necessarily the same as defined by
    // primitiveMesh::cellShapes().
    // Thus do the lookup for registered primitive types.
    // Finally, remap the cellModel face number to the ProstarFaceId
    if (prostarShapeLookup_.found(mapIndex))
    {
        const faceList sFaces = shape.faces();
        forAll(sFaces, sFaceI)
        {
            if (faces[faceI] == sFaces[sFaceI])
            {
                cellFaceId = sFaceI;

                if (ccm::debug > 1)
                {
                    Info << " FoamCellFace: " << sFaceI;
                }

                break;
            }
        }

        mapIndex = prostarShapeLookup_[mapIndex];
        cellFaceId = STARCDCore::foamToStarFaceAddr[mapIndex][cellFaceId];

        if (ccm::debug > 1)
        {
            Info<< " ProstarShape: " << mapIndex
                << " ProstarFaceId: " << cellFaceId + 1
                << endl;
        }
    }

    // PROSTAR used 1-based indices
    return cellFaceId + 1;
}


// write the faces, the number of vertices appears before each entry
void Foam::ccm::writer::writeFaces
(
    const ccmID& nodeId,
    const ccmID& mapId,
    bool  isBoundary,
    label size,
    label start
) const
{
    if (globalState_->hasError() || size <= 0)
    {
        return;
    }

    CCMIOEntity nodeType;
    if (isBoundary)
    {
        nodeType = kCCMIOBoundaryFaces;
    }
    else
    {
        nodeType = kCCMIOInternalFaces;
    }

    const faceList&  faces = mesh_.faces();
    const labelList& owner = mesh_.faceOwner();
    const labelList& neigh = mesh_.faceNeighbour();

    // 1. write vertices to define faces

    // Determine allocation size
    label streamSize = size;
    for (label faceI = 0; faceI < size; ++faceI)
    {
        streamSize += faces[faceI + start].size();
    }

    // Build a vertex stream
    List<int> ccmStream(streamSize);
    streamSize = 0;
    for (label faceI = 0; faceI < size; ++faceI)
    {
        const labelList& vrtList = faces[faceI + start];

        ccmStream[streamSize++] = vrtList.size();
        forAll(vrtList, i)
        {
            ccmStream[streamSize++] = vrtList[i] + 1;
        }
    }

    if (ccm::debug)
    {
        Info<<"CCMIOWriteFaces()  size:" << size << " streamSize:"
            << streamSize << endl;
    }

    CCMIOWriteFaces
    (
        &(globalState_->error),
        nodeId,
        nodeType,
        mapId,
        streamSize,
        ccmStream.cdata(),
        kCCMIOStart,
        kCCMIOEnd
    );

    if (globalState_->hasError())
    {
        return;
    }

    // 2. write face owner/neighbours
    if (isBoundary)
    {
        streamSize = size;
        if (ccmStream.size() < streamSize)
        {
            ccmStream.setSize(streamSize);
        }

        for (int faceI = 0; faceI < size; ++faceI)
        {
            ccmStream[faceI] = owner[faceI + start] + 1;
        }
    }
    else
    {
        // kCCMIOInternalFaces
        streamSize = 2 * size;
        if (ccmStream.size() < streamSize)
        {
            ccmStream.setSize(streamSize);
        }

        for (int faceI = 0; faceI < size; ++faceI)
        {
            ccmStream[faceI*2]   = owner[faceI + start] + 1;
            ccmStream[faceI*2+1] = neigh[faceI + start] + 1;
        }
    }

    if (ccm::debug)
    {
        Info<<"CCMIOWriteFaceCells()  size:" << size;
        if (isBoundary)
        {
            Info<< " boundary";
        }
        else
        {
            Info<<" internal";
        }
        Info<<"Faces"<<endl;
    }

    CCMIOWriteFaceCells
    (
        &(globalState_->error),
        nodeId,
        nodeType,
        mapId,
        ccmStream.cdata(),
        kCCMIOStart,
        kCCMIOEnd
    );

    // determine corresponding ProstarFaceId
    if (isBoundary)
    {
        streamSize = size;
        if (ccmStream.size() < streamSize)
        {
            ccmStream.setSize(streamSize);
        }

        for (int faceIdx = 0; faceIdx < size; ++faceIdx)
        {
            label faceI = faceIdx + start;
            // boundary face - owner only
            ccmStream[faceIdx] = prostarCellFaceId(owner[faceI], faceI);
        }

        CCMIOWriteOpt1i
        (
            &(globalState_->error),
            nodeId,
            "ProstarFaceId",
            size,
            ccmStream.cdata(),
            kCCMIOStart,
            kCCMIOEnd
        );
    }
    else
    {
        // kCCMIOInternalFaces
        streamSize = 2 * size;
        if (ccmStream.size() < streamSize)
        {
            ccmStream.setSize(streamSize);
        }

        for (int faceIdx = 0; faceIdx < size; ++faceIdx)
        {
            label faceI = faceIdx + start;
            ccmStream[faceIdx*2]   = prostarCellFaceId(owner[faceI], faceI);
            ccmStream[faceIdx*2+1] = prostarCellFaceId(neigh[faceI], faceI);
        }

        CCMIOWriteOpt2i
        (
            &(globalState_->error),
            nodeId,
            "ProstarFaceId",
            size,
            2,
            ccmStream.cdata(),
            kCCMIOStart,
            kCCMIOEnd
        );
    }
}


// writeVertices
// 1) write the vertex map (starting with 1)
// 2) write the vertex data
//
void Foam::ccm::writer::writeVertices
(
    const ccmID& verticesNode
) const
{
    const pointField& pts = mesh_.points();

    Info<< "writing points: " << pts.size() << endl;

    // 1. mapping data array index to the vertex Id (starting with 1)
    ccmID vertexMap;
    addLinearMap
    (
        "vertex Map",
        vertexMap,
        pts.size()
    );

    // 2. write vertices - scale [m] -> [mm] for consistency with PROSTAR
    //    but is probably immaterial
    const float scaling(0.001);   // to recover meters
    List<float> vrts(3*pts.size());
    forAll(pts, i)
    {
        vrts[3*i]    = 1000.0 * pts[i].x();
        vrts[3*i+1]  = 1000.0 * pts[i].y();
        vrts[3*i+2]  = 1000.0 * pts[i].z();
    }

    CCMIOWriteVerticesf
    (
        &(globalState_->error),
        verticesNode,
        3, scaling,
        vertexMap,
        vrts.cdata(),
        kCCMIOStart,
        kCCMIOEnd
    );
    assertNoError("writing 'Vertices' node");
}


// writeInternalFaces
// 1) write the face map (starting with 1)
// 2) write owner/neighbour
//
void Foam::ccm::writer::writeInternalFaces
(
    const ccmID& topoNode
) const
{
    label nFaces = mesh_.nInternalFaces();

    Info<< "writing internalFaces: " << nFaces << endl;
    if (nFaces > 0)
    {
        ccmID nodeId;
        CCMIONewEntity
        (
            &(globalState_->error),
            topoNode,
            kCCMIOInternalFaces,
            nullptr,
            &nodeId
        );
        assertNoError("creating internalFaces node");

        writeFaces
        (
            nodeId,
            maps_->internalFaces,
            false,
            nFaces
        );
        assertNoError("writing internalFaces");
    }
}


// writeBoundaryFaces:
// - write faces with owner
void Foam::ccm::writer::writeBoundaryFaces
(
    const ccmID& topoNode
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // 4. write boundary faces
    Info<< "writing boundaryFaces:" << flush;

    label defaultId = findDefaultBoundary();

    // write Default_Boundary_Region as BoundaryFaces-0
    if (defaultId >= 0)
    {
        Info<< " " << patches[defaultId].size() << flush;

        ccmID nodeId;
        CCMIONewIndexedEntity
        (
            &(globalState_->error),
            topoNode,
            kCCMIOBoundaryFaces,
            0,
            nullptr,
            &nodeId
        );
        assertNoError
        (
            "creating boundaryFaces node patch: " + patches[defaultId].name()
        );

        writeFaces
        (
            nodeId,
            maps_->boundary[defaultId],
            true,
            patches[defaultId].size(),
            patches[defaultId].start()
        );
        assertNoError
        (
            "writing boundaryFaces patch: " + patches[defaultId].name()
        );
    }

    //
    // write boundary faces - skip Default_Boundary_Region
    //
    forAll(patches, patchI)
    {
        label regionId = patchI;
        if (regionId == defaultId)
        {
            continue;  // skip - already written
        }
        else if (defaultId == -1 || regionId < defaultId)
        {
            regionId++;
        }

        Info<< " " << patches[patchI].size() << flush;

        if (patches[patchI].size() > 0)
        {
            ccmID nodeId;
            CCMIONewIndexedEntity
            (
                &(globalState_->error),
                topoNode,
                kCCMIOBoundaryFaces,
                regionId,
                nullptr,
                &nodeId
            );
            assertNoError
            (
                "creating boundaryFaces node patch: " + patches[patchI].name()
            );

            writeFaces
            (
                nodeId,
                maps_->boundary[patchI],
                true,
                patches[patchI].size(),
                patches[patchI].start()
            );
            assertNoError
            (
                "writing boundaryFaces patch: " + patches[patchI].name()
            );
        }
    }

    Info<< endl;
}


// writeCells:
// - write faces with owner/neighbour
// - write interfaces
void Foam::ccm::writer::writeCells
(
    const ccmID& topoNode
)
{
    Info<< "writing cells: " << mesh_.nCells() << endl;

    // 1. cellTableId
    //    - if possible, read from constant/polyMesh/cellTableId
    //    - otherwise use cellZone information
    List<int> mapData(mesh_.nCells(), -1);
    bool useCellZones = false;

    IOList<label> ioList
    (
        IOobject
        (
            "cellTableId",
            "constant",
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );


    if (ioList.headerOk())
    {
        if (ioList.size() == mesh_.nCells())
        {
            mapData.transfer(ioList);
        }
        else
        {
            WarningInFunction
                << ioList.objectPath() << endl
                << "    Has incorrect number of cells: "
                << ioList.size() << " instead of " << mesh_.nCells()
                << " -    use cellZone information instead"
                << endl;

            ioList.clear();
            useCellZones = true;
        }
    }
    else
    {
        useCellZones = true;
    }


    if (useCellZones)
    {
        if (!cellTable_.size())
        {
            Info<< "created cellTable from cellZones" << endl;
            cellTable_ = mesh_;
        }

        // track if there are unzoned cells
        label nUnzoned = mesh_.nCells();

        // get the cellZone <-> cellTable correspondence
        Info<< "matching cellZones to cellTable" << endl;

        forAll(mesh_.cellZones(), zoneI)
        {
            const cellZone& cZone = mesh_.cellZones()[zoneI];
            if (cZone.size())
            {
                nUnzoned -= cZone.size();

                label tableId = cellTable_.findIndex(cZone.name());
                if (tableId < 0)
                {
                    dictionary dict;

                    dict.add("Label", cZone.name());
                    dict.add("MaterialType", "fluid");
                    tableId = cellTable_.append(dict);
                }

                for (auto id : cZone)
                {
                    mapData[id] = tableId;
                }
            }
        }

        if (nUnzoned)
        {
            dictionary dict;

            dict.add("Label", "__unzonedCells__");
            dict.add("MaterialType", "fluid");
            label tableId = cellTable_.append(dict);

            for (auto& id : mapData)
            {
                if (id < 0)
                {
                    id = tableId;
                }
            }
        }
    }

    // 2. mapping data array index to the cell Id (starting with 1)
    ccmID cellsNode;
    CCMIONewEntity
    (
        &(globalState_->error),
        topoNode,
        kCCMIOCells,
        nullptr,
        &cellsNode
    );
    assertNoError("creating 'Cells' node");

    CCMIOWriteCells
    (
        &(globalState_->error),
        cellsNode,
        maps_->cells,
        mapData.data(),
        kCCMIOStart, kCCMIOEnd
    );
    assertNoError("writing 'Cells' node");

    // 3. add cell topology information, if possible
    const cellShapeList& shapes = mesh_.cellShapes();
    forAll(shapes, cellI)
    {
        label mapIndex = shapes[cellI].model().index();

        // A registered primitive type
        if (prostarShapeLookup_.found(mapIndex))
        {
            mapData[cellI] = prostarShapeLookup_[mapIndex];
        }
        else
        {
            mapData[cellI] = STARCDCore::starcdPoly;  // Treat as polyhedral
        }
    }

    CCMIOWriteOpt1i
    (
        &(globalState_->error),
        cellsNode,
        "CellTopologyType",
        mesh_.nCells(),
        mapData.cdata(),
        kCCMIOStart, kCCMIOEnd
    );

    // 4. write interfaces
    writeInterfaces(cellsNode);
}


// write interfaces
//   1) PROSTAR baffles
//
void Foam::ccm::writer::writeInterfaces
(
    const ccmID& cellsNode
) const
{
    IOList<labelList> interfaces
    (
        IOobject
        (
            "interfaces",
            "constant",
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    if (interfaces.headerOk() && interfaces.size() > 0)
    {
        List<int> mapData(2 * interfaces.size());

        forAll(interfaces, i)
        {
            mapData[i*2]   = interfaces[i][0];
            mapData[i*2+1] = interfaces[i][1];
        }

        ccmID nodeId;
        if
        (
            CCMIONewEntity
            (
                &(globalState_->error),
                cellsNode,
                kCCMIOInterfaces,
                0,
                &nodeId
            )
         != kCCMIONoErr

         || CCMIOWriteOpt2i
            (
                &(globalState_->error),
                nodeId,
                "FaceIds",
                interfaces.size(),
                2,
                mapData.cdata(), kCCMIOStart, kCCMIOEnd
            )
         != kCCMIONoErr
        )
        {
            assertNoError("writing interfaces 'FaceIds'");
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ccm::writer::writeGeometry()
{
    // use first ("default") state node
    ccmID stateNode;
    // int stateI = 0;

    // create "default" State
    CCMIONewState
    (
        &(globalState_->error),
        (globalState_->root),
        "default",
        nullptr,
        nullptr,
        &stateNode
    );
    assertNoError("could not create default state");

    // use first processor - create a new one
    ccmID processorNode;
    int procI = 0;

    if
    (
        CCMIONextEntity
        (
            nullptr,
            stateNode,
            kCCMIOProcessor,
            &procI,
            &processorNode
        )
     != kCCMIONoErr
    )
    {
        CCMIONewEntity
        (
            &(globalState_->error),
            stateNode,
            kCCMIOProcessor,
            nullptr,
            &processorNode
        );
        assertNoError("could not create processor node");
    }

    // remove old data (if it exists)
    CCMIOClearProcessor
    (
        nullptr, stateNode, processorNode,
        true,   // Clear vertices
        true,   // Clear topology
        true,   // Clear initial field
        true,   // Clear solution
        true    // Clear lagrangian
    );

#if 0
    CCMIOWriteOptstr
    (
        nullptr,
        processorNode,
        "CreatingProgram",
        "ccm::writer"
    );
#endif

    //
    // create vertices and topology nodes
    //
    ccmID verticesNode, topoNode;
    if
    (
        CCMIONewEntity
        (
            &(globalState_->error),
            (globalState_->root),
            kCCMIOVertices,
            nullptr,
            &verticesNode
        )
     == kCCMIONoErr

     && CCMIONewEntity
        (
            &(globalState_->error),
            (globalState_->root),
            kCCMIOTopology,
            nullptr,
            &topoNode
        )
     == kCCMIONoErr
    )
    {
        writeVertices(verticesNode);
        writeInternalFaces(topoNode);
        writeBoundaryFaces(topoNode);
        writeCells(topoNode);
        writeProblem(stateNode);
    }

    // Now we have the mesh (vertices and topology),
    // we can write out the processor information.
    CCMIOWriteProcessor
    (
        &(globalState_->error),
        processorNode,
        nullptr, &verticesNode, // no verticesFile, write verticesNode
        nullptr, &topoNode,     // no topologyFile, write topoNode
        nullptr, nullptr,       // initialField unchanged
        nullptr, nullptr        // no solutionFile, solution unchanged
    );
    assertNoError("Error after writing geometry processor");

}


// ************************************************************************* //
