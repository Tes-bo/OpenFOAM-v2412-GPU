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

#include "PatchTools.H"
#include "polyMesh.H"
#include "indirectPrimitivePatch.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::tmp<Foam::pointField>
Foam::PatchTools::pointNormals
(
    const polyMesh& mesh,
    const PrimitivePatch<FaceList, PointField>& p,
    const bitSet& pFlip
)
{
    const globalMeshData& globalData = mesh.globalData();
    const indirectPrimitivePatch& coupledPatch = globalData.coupledPatch();
    const Map<label>& coupledPatchMP = coupledPatch.meshPointMap();
    const mapDistribute& map = globalData.globalPointSlavesMap();
    const globalIndexAndTransform& transforms =
        globalData.globalTransforms();


    // Combine normals. Note: do on all master points. Cannot just use
    // patch points since the master point does not have to be on the
    // patch!

    pointField coupledPointNormals(map.constructSize(), Zero);

    {
        // Collect local pointFaces (sized on patch points only)
        List<List<point>> pointFaceNormals(map.constructSize());
        forAll(p.meshPoints(), patchPointi)
        {
            const label meshPointi = p.meshPoints()[patchPointi];

            const auto fnd = coupledPatchMP.cfind(meshPointi);
            if (fnd.good())
            {
                const label coupledPointi = fnd.val();

                List<point>& pNormals = pointFaceNormals[coupledPointi];
                const labelList& pFaces = p.pointFaces()[patchPointi];
                pNormals.setSize(pFaces.size());
                forAll(pFaces, i)
                {
                    const label facei = pFaces[i];
                    const vector& n = p.faceNormals()[facei];
                    pNormals[i] = ((pFlip.empty() || !pFlip[facei]) ? n : -n);
                }
            }
        }


        // Pull remote data into local slots
        map.distribute
        (
            transforms,
            pointFaceNormals,
            mapDistribute::transform()
        );


        // Combine all face normals (-local, -remote,untransformed,
        //  -remote,transformed)

        const labelListList& slaves = globalData.globalPointSlaves();
        const labelListList& transformedSlaves =
            globalData.globalPointTransformedSlaves();

        forAll(slaves, coupledPointi)
        {
            const labelList& slaveSlots = slaves[coupledPointi];
            const labelList& transformedSlaveSlots =
                transformedSlaves[coupledPointi];

            point& n = coupledPointNormals[coupledPointi];

            // Local entries
            const List<point>& local = pointFaceNormals[coupledPointi];

            label nFaces =
                local.size()
              + slaveSlots.size()
              + transformedSlaveSlots.size();

            n = sum(local);

            // Add any remote face normals
            forAll(slaveSlots, i)
            {
                n += sum(pointFaceNormals[slaveSlots[i]]);
            }
            forAll(transformedSlaveSlots, i)
            {
                n += sum(pointFaceNormals[transformedSlaveSlots[i]]);
            }

            if (nFaces >= 1)
            {
                n /= mag(n)+VSMALL;
            }

            // Put back into slave slots
            forAll(slaveSlots, i)
            {
                coupledPointNormals[slaveSlots[i]] = n;
            }
            forAll(transformedSlaveSlots, i)
            {
                coupledPointNormals[transformedSlaveSlots[i]] = n;
            }
        }


        // Send back
        map.reverseDistribute
        (
            transforms,
            coupledPointNormals.size(),
            coupledPointNormals,
            mapDistribute::transform()
        );
    }


    // 1. Start off with local normals (note:without calculating pointNormals
    //    to avoid them being stored)

    auto textrudeN = tmp<pointField>::New(p.nPoints(), Zero);
    auto& extrudeN = textrudeN.ref();
    {
        const faceList& localFaces = p.localFaces();
        const vectorField& faceNormals = p.faceNormals();

        forAll(localFaces, facei)
        {
            const face& f = localFaces[facei];
            const vector& n = faceNormals[facei];
            forAll(f, fp)
            {
                extrudeN[f[fp]] += ((pFlip.empty() || !pFlip[facei]) ? n : -n);
            }
        }
        extrudeN /= mag(extrudeN)+VSMALL;
    }


    // 2. Override patch normals on coupled points
    forAll(p.meshPoints(), patchPointi)
    {
        const label meshPointi = p.meshPoints()[patchPointi];

        const auto fnd = coupledPatchMP.cfind(meshPointi);
        if (fnd.good())
        {
            const label coupledPointi = fnd.val();
            extrudeN[patchPointi] = coupledPointNormals[coupledPointi];
        }
    }

    return textrudeN;
}


template<class FaceList, class PointField>
Foam::tmp<Foam::pointField>
Foam::PatchTools::edgeNormals
(
    const polyMesh& mesh,
    const PrimitivePatch<FaceList, PointField>& p,
    const labelList& patchEdges,
    const labelList& coupledEdges,
    const bitSet& pFlip
)
{
    // 1. Start off with local normals

    auto tedgeNormals = tmp<pointField>::New(p.nEdges(), Zero);
    auto& edgeNormals = tedgeNormals.ref();

    {
        const labelListList& edgeFaces = p.edgeFaces();
        const vectorField& faceNormals = p.faceNormals();

        forAll(edgeFaces, edgei)
        {
            const labelList& eFaces = edgeFaces[edgei];
            for (const label facei : eFaces)
            {
                const vector& n = faceNormals[facei];
                edgeNormals[edgei] +=
                (
                    (pFlip.empty() || !pFlip[facei])
                  ? n
                  : -n
                );
            }
        }
        edgeNormals /= mag(edgeNormals)+VSMALL;
    }



    const globalMeshData& globalData = mesh.globalData();
    const mapDistribute& map = globalData.globalEdgeSlavesMap();


    // Convert patch-edge data into cpp-edge data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Construct with all data in consistent orientation
    pointField cppEdgeData(map.constructSize(), Zero);

    forAll(patchEdges, i)
    {
        label patchEdgeI = patchEdges[i];
        label coupledEdgeI = coupledEdges[i];
        cppEdgeData[coupledEdgeI] = edgeNormals[patchEdgeI];
    }


    // Synchronise
    // ~~~~~~~~~~~

    globalData.syncData
    (
        cppEdgeData,
        globalData.globalEdgeSlaves(),
        globalData.globalEdgeTransformedSlaves(),
        map,
        globalData.globalTransforms(),
        plusEqOp<point>(),              // add since normalised later on
        mapDistribute::transform()
    );
    cppEdgeData /= mag(cppEdgeData)+VSMALL;


    // Back from cpp-edge to patch-edge data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(patchEdges, i)
    {
        label patchEdgeI = patchEdges[i];
        label coupledEdgeI = coupledEdges[i];
        edgeNormals[patchEdgeI] = cppEdgeData[coupledEdgeI];
    }

    return tedgeNormals;
}


// ************************************************************************* //
