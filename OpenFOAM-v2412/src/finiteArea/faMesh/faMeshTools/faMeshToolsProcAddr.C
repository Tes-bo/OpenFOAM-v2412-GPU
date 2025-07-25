/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "faMeshTools.H"
#include "BitOps.H"
#include "fileOperation.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "IOmapDistributePolyMesh.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Create a reconstruct map.

static mapDistributePolyMesh createReconstructMap
(
    const faMesh& mesh,
    const autoPtr<faMesh>& baseMeshPtr,
    const labelUList& faceProcAddr,
    const labelUList& edgeProcAddr,
    const labelUList& pointProcAddr,
    const labelUList& boundaryProcAddr
)
{
    const label nOldPoints = mesh.nPoints();
    const label nOldFaces = mesh.nFaces();
    const label nOldEdges = mesh.nEdges();

    ///Pout<< "old sizes"
    ///    << " points:" << nOldPoints
    ///    << " faces:" << nOldFaces
    ///    << " edges:" << nOldEdges << nl;

    const faBoundaryMesh& oldBndMesh = mesh.boundary();
    labelList oldPatchStarts(oldBndMesh.patchStarts());

    // Patches: purge -1 entries
    labelList patchProcAddr
    (
        IndirectList<label>::subset_if
        (
            boundaryProcAddr,
            labelRange::ge0()
        )
    );


    labelListList faceSubMap(Pstream::nProcs());
    faceSubMap[Pstream::masterNo()] = identity(nOldFaces);

    labelListList edgeSubMap(Pstream::nProcs());
    edgeSubMap[Pstream::masterNo()] = identity(nOldEdges);

    labelListList pointSubMap(Pstream::nProcs());
    pointSubMap[Pstream::masterNo()] = identity(nOldPoints);

    labelListList patchSubMap(Pstream::nProcs());
    patchSubMap[Pstream::masterNo()] = patchProcAddr;


    // Gather addressing on the master
    labelListList faceAddressing(Pstream::nProcs());
    faceAddressing[Pstream::myProcNo()] = faceProcAddr;
    Pstream::gatherList(faceAddressing);

    labelListList edgeAddressing(Pstream::nProcs());
    edgeAddressing[Pstream::myProcNo()] = edgeProcAddr;
    Pstream::gatherList(edgeAddressing);

    labelListList pointAddressing(Pstream::nProcs());
    pointAddressing[Pstream::myProcNo()] = pointProcAddr;
    Pstream::gatherList(pointAddressing);

    labelListList patchAddressing(Pstream::nProcs());
    patchAddressing[Pstream::myProcNo()] = patchProcAddr;
    Pstream::gatherList(patchAddressing);


    // NB: can only have a reconstruct on master!
    if (Pstream::master() && baseMeshPtr && baseMeshPtr->nFaces())
    {
        const faMesh& baseMesh = *baseMeshPtr;

        const label nNewPoints = baseMesh.nPoints();
        const label nNewFaces = baseMesh.nFaces();
        const label nNewEdges = baseMesh.nEdges();
        const label nNewPatches = baseMesh.boundary().size();

        /// Pout<< "new sizes"
        ///     << " points:" << nNewPoints
        ///     << " faces:" << nNewFaces
        ///     << " edges:" << nNewEdges << nl;

        mapDistribute faFaceMap
        (
            nNewFaces,
            std::move(faceSubMap),
            std::move(faceAddressing),
            false,  // subHasFlip
            false   // constructHasFlip
        );

        mapDistribute faEdgeMap
        (
            nNewEdges,
            std::move(edgeSubMap),
            std::move(edgeAddressing),
            false,  // subHasFlip
            false   // constructHasFlip
        );

        mapDistribute faPointMap
        (
            nNewPoints,
            std::move(pointSubMap),
            std::move(pointAddressing)
        );

        mapDistribute faPatchMap
        (
            nNewPatches,
            std::move(patchSubMap),
            std::move(patchAddressing)
        );

        return mapDistributePolyMesh
        (
            // Mesh before changes
            nOldPoints,
            nOldEdges,          // area: nOldEdges (volume: nOldFaces)
            nOldFaces,          // area: nOldFaces (volume: nOldCells)

            std::move(oldPatchStarts),
            labelList(),        // oldPatchNMeshPoints [unused]

            mapDistribute(std::move(faPointMap)),
            mapDistribute(std::move(faEdgeMap)), // edgeMap (volume: faceMap)
            mapDistribute(std::move(faFaceMap)), // faceMap (volume: cellMap)
            mapDistribute(std::move(faPatchMap))
        );
    }
    else
    {
        mapDistribute faFaceMap
        (
            0,  // nNewFaces
            std::move(faceSubMap),
            labelListList(Pstream::nProcs()),   // constructMap
            false,  // subHasFlip
            false   // constructHasFlip
        );

        mapDistribute faEdgeMap
        (
            0,  // nNewEdges
            std::move(edgeSubMap),
            labelListList(Pstream::nProcs()),   // constructMap
            false,  // subHasFlip
            false   // constructHasFlip
        );

        mapDistribute faPointMap
        (
            0,  // nNewPoints
            std::move(pointSubMap),
            labelListList(Pstream::nProcs())    // constructMap
        );

        mapDistribute faPatchMap
        (
            0,  // nNewPatches
            std::move(patchSubMap),
            labelListList(Pstream::nProcs())    // constructMap
        );

        return mapDistributePolyMesh
        (
            // Mesh before changes
            nOldPoints,
            nOldEdges,          // area: nOldEdges (volume: nOldFaces)
            nOldFaces,          // area: nOldFaces (volume: nOldCells)

            std::move(oldPatchStarts),
            labelList(),        // oldPatchNMeshPoints [unused]

            mapDistribute(std::move(faPointMap)),
            mapDistribute(std::move(faEdgeMap)), // edgeMap (volume: faceMap)
            mapDistribute(std::move(faFaceMap)), // faceMap (volume: cellMap)
            mapDistribute(std::move(faPatchMap))
        );
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::mapDistributePolyMesh
Foam::faMeshTools::readProcAddressing
(
    const faMesh& mesh,
    const autoPtr<faMesh>& baseMeshPtr
)
{
    IOobject ioAddr
    (
        "procAddressing",
        mesh.facesInstance(),
        faMesh::meshSubDir,
        mesh.thisDb(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    //if (ioAddr.typeHeaderOk<labelIOList>(true))
    //{
    //    Pout<< "Reading addressing from " << io.name() << " at "
    //        << mesh.facesInstance() << nl << endl;
    //    distMap.reset(new IOmapDistributePolyMesh(io));
    //}
    //else

    {
        Info<< "Reading (face|edge|face|point|boundary)ProcAddressing from "
            << mesh.facesInstance().c_str() << '/'
            << faMesh::meshSubDir << nl << endl;

        ioAddr.rename("faceProcAddressing");
        labelIOList faceProcAddressing(ioAddr, Zero);

        ioAddr.rename("edgeProcAddressing");
        labelIOList edgeProcAddressing(ioAddr, Zero);

        ioAddr.rename("pointProcAddressing");
        labelIOList pointProcAddressing(ioAddr, Zero);

        ioAddr.rename("boundaryProcAddressing");
        labelIOList boundaryProcAddressing(ioAddr, Zero);

        if
        (
            mesh.nFaces() != faceProcAddressing.size()
         || mesh.nEdges() != edgeProcAddressing.size()
         || mesh.nPoints() != pointProcAddressing.size()
         || mesh.boundary().size() != boundaryProcAddressing.size()
        )
        {
            FatalErrorInFunction
                << "Read addressing inconsistent with mesh sizes" << nl
                << "faces:" << mesh.nFaces()
                << " addressing:" << faceProcAddressing.objectRelPath()
                << " size:" << faceProcAddressing.size() << nl
                << "edges:" << mesh.nEdges()
                << " addressing:" << edgeProcAddressing.objectRelPath()
                << " size:" << edgeProcAddressing.size() << nl
                << "points:" << mesh.nPoints()
                << " addressing:" << pointProcAddressing.objectRelPath()
                << " size:" << pointProcAddressing.size()
                << "patches:" << mesh.boundary().size()
                << " addressing:" << boundaryProcAddressing.objectRelPath()
                << " size:" << boundaryProcAddressing.size()
                << exit(FatalError);
        }

        return createReconstructMap
        (
            mesh,
            baseMeshPtr,
            faceProcAddressing,
            edgeProcAddressing,
            pointProcAddressing,
            boundaryProcAddressing
        );
    }
}


void Foam::faMeshTools::writeProcAddressing
(
    const faMesh& mesh,
    const mapDistributePolyMesh& map,
    const bool decompose,
    refPtr<fileOperation>& writeHandler,
    const faMesh* procMesh
)
{
    Info<< "Writing ("
        << (decompose ? "decompose" : "reconstruct")
        << ") procAddressing files to "
        << mesh.facesInstance().c_str() << '/'
        << faMesh::meshSubDir << endl;

    IOobject ioAddr
    (
        "procAddressing",
        mesh.facesInstance(),
        faMesh::meshSubDir,
        (procMesh && !decompose ? procMesh->thisDb() : mesh.thisDb()),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );


    // faceProcAddressing (faMesh)
    ioAddr.rename("faceProcAddressing");
    labelIOList faceMap(ioAddr, Zero);

    // edgeProcAddressing (faMesh)
    ioAddr.rename("edgeProcAddressing");
    labelIOList edgeMap(ioAddr, Zero);

    // pointProcAddressing (faMesh)
    ioAddr.rename("pointProcAddressing");
    labelIOList pointMap(ioAddr, Zero);

    // boundaryProcAddressing (faMesh)
    ioAddr.rename("boundaryProcAddressing");
    labelIOList patchMap(ioAddr, Zero);

    if (decompose)
    {
        // Decompose
        // - forward map:  [undecomposed] -> [decomposed]

        // area:faces (volume:cells)
        faceMap = identity(map.nOldCells());
        map.cellMap().distribute(faceMap);

        // area:edges (volume:faces)
        edgeMap = identity(map.nOldFaces());
        map.faceMap().distribute(edgeMap);

        pointMap = identity(map.nOldPoints());
        map.distributePointData(pointMap);

        patchMap = identity(map.patchMap().constructSize());
        map.patchMap().mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            label(-1),  // nullValue for new patches...
            patchMap,
            flipOp()    // negate op
        );
    }
    else    // reconstruct
    {
        // Reconstruct
        // - reverse map:  [undecomposed] <- [decomposed]

        // area:faces (volume:cells)
        faceMap = identity(mesh.nFaces());
        map.cellMap().reverseDistribute(map.nOldCells(), faceMap);

        // area:edges (volume:faces)
        edgeMap = identity(mesh.patch().nEdges());
        map.faceMap().reverseDistribute(map.nOldFaces(), edgeMap);

        pointMap = identity(mesh.nPoints());
        map.pointMap().reverseDistribute(map.nOldPoints(), pointMap);

        patchMap = identity(mesh.boundary().size());
        map.patchMap().mapDistributeBase::reverseDistribute
        (
            Pstream::commsTypes::nonBlocking,
            map.oldPatchSizes().size(),
            label(-1),  // nullValue for unmapped patches...
            patchMap
        );
    }

    auto oldHandler = fileOperation::fileHandler(writeHandler);


    // If we want procAddressing, need to manually write it ourselves
    // since it was not registered anywhere

    IOmapDistributePolyMeshRef procAddrMap
    (
        IOobject
        (
            "procAddressing",
            mesh.facesInstance(),
            faMesh::meshSubDir,
            mesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        map
    );


    if (decompose)
    {
        // Write into proc directories
        procAddrMap.write();
    }
    else
    {
        // Reconstruct: "procAddressing" only meaningful for rank 0
        // and written into base (serial) location (if at all).

        if (UPstream::master())
        {
            const bool oldParRun = UPstream::parRun(false);
            procAddrMap.write();
            UPstream::parRun(oldParRun);
        }
    }


    const bool faceOk = faceMap.write();
    const bool edgeOk = edgeMap.write();
    const bool pointOk = pointMap.write();
    const bool patchOk = patchMap.write();

    writeHandler = fileOperation::fileHandler(oldHandler);

    if (!edgeOk || !faceOk || !pointOk || !patchOk)
    {
        WarningInFunction
            << "Failed to write some of "
            << faceMap.objectRelPath() << ", "
            << edgeMap.objectRelPath() << ", "
            << pointMap.objectRelPath() << ", "
            << patchMap.objectRelPath() << endl;
    }
}


// ************************************************************************* //
