/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Application
    syncToolsTest

Description
    Test some functionality in syncTools.

\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "Random.H"
#include "PackedList.H"
#include "flipOp.H"
#include "pointList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void testPackedList(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing PackedList synchronisation." << endl;

    {
        PackedList<3> bits(mesh.nEdges());
        forAll(bits, i)
        {
            bits.set(i, rndGen.position<label>(0,3));
        }

        labelList edgeValues(mesh.nEdges());
        forAll(bits, i)
        {
            edgeValues[i] = bits.get(i);
        }

        PackedList<3> maxBits(bits);
        labelList maxEdgeValues(edgeValues);

        syncTools::syncEdgeList(mesh, bits, minEqOp<unsigned int>(), 0);
        syncTools::syncEdgeList(mesh, edgeValues, minEqOp<label>(), 0);

        syncTools::syncEdgeList(mesh, maxBits, maxEqOp<unsigned int>(), 0);
        syncTools::syncEdgeList
        (
            mesh,
            maxEdgeValues,
            maxEqOp<label>(),
            0
        );

        forAll(bits, i)
        {
            if
            (
                edgeValues[i] != label(bits.get(i))
             || maxEdgeValues[i] != label(maxBits.get(i))
            )
            {
                FatalErrorInFunction
                    << "edge:" << i
                    << " minlabel:" << edgeValues[i]
                    << " minbits:" << bits.get(i)
                    << " maxLabel:" << maxEdgeValues[i]
                    << " maxBits:" << maxBits.get(i)
                    << exit(FatalError);
            }
        }
    }

    {
        PackedList<3> bits(mesh.nPoints());
        forAll(bits, i)
        {
            bits.set(i, rndGen.position<label>(0,3));
        }

        labelList pointValues(mesh.nPoints());
        forAll(bits, i)
        {
            pointValues[i] = bits.get(i);
        }

        PackedList<3> maxBits(bits);
        labelList maxPointValues(pointValues);

        syncTools::syncPointList(mesh, bits, minEqOp<unsigned int>(), 0);
        syncTools::syncPointList(mesh, pointValues, minEqOp<label>(), 0);

        syncTools::syncPointList(mesh, maxBits, maxEqOp<unsigned int>(), 0);
        syncTools::syncPointList
        (
            mesh,
            maxPointValues,
            maxEqOp<label>(),
            0
        );

        forAll(bits, i)
        {
            if
            (
                pointValues[i] != label(bits.get(i))
             || maxPointValues[i] != label(maxBits.get(i))
            )
            {
                FatalErrorInFunction
                    << "point:" << i
                    << " at:" << mesh.points()[i]
                    << " minlabel:" << pointValues[i]
                    << " minbits:" << bits.get(i)
                    << " maxLabel:" << maxPointValues[i]
                    << " maxBits:" << maxBits.get(i)
                    << exit(FatalError);
            }
        }
    }

    {
        PackedList<3> bits(mesh.nFaces());
        forAll(bits, facei)
        {
            bits.set(facei, rndGen.position<label>(0,3));
        }

        labelList faceValues(mesh.nFaces());
        forAll(bits, facei)
        {
            faceValues[facei] = bits.get(facei);
        }

        PackedList<3> maxBits(bits);
        labelList maxFaceValues(faceValues);

        syncTools::syncFaceList(mesh, bits, minEqOp<unsigned int>());
        syncTools::syncFaceList(mesh, faceValues, minEqOp<label>());

        syncTools::syncFaceList(mesh, maxBits, maxEqOp<unsigned int>());
        syncTools::syncFaceList(mesh, maxFaceValues, maxEqOp<label>());

        forAll(bits, facei)
        {
            if
            (
                faceValues[facei] != label(bits.get(facei))
             || maxFaceValues[facei] != label(maxBits.get(facei))
            )
            {
                FatalErrorInFunction
                    << "face:" << facei
                    << " minlabel:" << faceValues[facei]
                    << " minbits:" << bits.get(facei)
                    << " maxLabel:" << maxFaceValues[facei]
                    << " maxBits:" << maxBits.get(facei)
                    << exit(FatalError);
            }
        }
    }
}


void testSparseData(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing Map synchronisation." << endl;

    WarningInFunction
        << "Position test of sparse data only correct for cases without cyclics"
        << " with shared points." << endl;

    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nBoundaryFaces(),
            mesh.nInternalFaces()
        ),
        mesh.points()
    );
    const pointField& localPoints = allBoundary.localPoints();


    // Point data
    // ~~~~~~~~~~

    {
        // Create some data. Use slightly perturbed positions.
        Map<point> sparseData;
        pointField fullData(mesh.nPoints(), point(GREAT, GREAT, GREAT));

        forAll(localPoints, i)
        {
            const point pt = localPoints[i] + 1e-4*rndGen.sample01<vector>();

            label meshPointi = allBoundary.meshPoints()[i];

            sparseData.insert(meshPointi, pt);
            fullData[meshPointi] = pt;
        }

        syncTools::syncPointMap
        (
            mesh,
            sparseData,
            minMagSqrEqOp<point>()
            // true                    // apply separation
        );
        syncTools::syncPointList
        (
            mesh,
            fullData,
            minMagSqrEqOp<point>(),
            point(GREAT, GREAT, GREAT)
            // true                    // apply separation
        );

        // Compare.
        // 1. Is all fullData also present in sparseData and same value
        forAll(fullData, meshPointi)
        {
            const point& fullPt = fullData[meshPointi];

            if (fullPt != point(GREAT, GREAT, GREAT))
            {
                const point& sparsePt = sparseData[meshPointi];

                if (fullPt != sparsePt)
                {
                    FatalErrorInFunction
                        << "point:" << meshPointi
                        << " full:" << fullPt
                        << " sparse:" << sparsePt
                        << exit(FatalError);
                }
            }
        }

        // 2. Does sparseData contain more?
        forAllConstIters(sparseData, iter)
        {
            const label meshPointi = iter.key();
            const point& sparsePt = iter.val();
            const point& fullPt = fullData[meshPointi];

            if (fullPt != sparsePt)
            {
                FatalErrorInFunction
                    << "point:" << meshPointi
                    << " full:" << fullPt
                    << " sparse:" << sparsePt
                    << exit(FatalError);
            }
        }
    }


    // Edge data
    // ~~~~~~~~~

    {
        // Create some data. Use slightly perturbed positions.
        EdgeMap<point> sparseData;
        pointField fullData(mesh.nEdges(), point(GREAT, GREAT, GREAT));

        const edgeList& edges = allBoundary.edges();
        const labelList meshEdges = allBoundary.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        );

        forAll(edges, i)
        {
            const edge& e = edges[i];

            const point pt =
                e.centre(localPoints) + 1e-4*rndGen.sample01<vector>();

            label meshEdgeI = meshEdges[i];

            sparseData.insert(mesh.edges()[meshEdgeI], pt);
            fullData[meshEdgeI] = pt;
        }

        syncTools::syncEdgeMap
        (
            mesh,
            sparseData,
            minMagSqrEqOp<point>()
        );
        syncTools::syncEdgeList
        (
            mesh,
            fullData,
            minMagSqrEqOp<point>(),
            point(GREAT, GREAT, GREAT)
        );

        // Compare.
        // 1. Is all fullData also present in sparseData and same value
        forAll(fullData, meshEdgeI)
        {
            const point& fullPt = fullData[meshEdgeI];

            if (fullPt != point(GREAT, GREAT, GREAT))
            {
                const point& sparsePt = sparseData[mesh.edges()[meshEdgeI]];

                if (fullPt != sparsePt)
                {
                    FatalErrorInFunction
                        << "edge:" << meshEdgeI
                        << " points:" << mesh.edges()[meshEdgeI]
                        << " full:" << fullPt
                        << " sparse:" << sparsePt
                        << exit(FatalError);
                }
            }
        }

        // 2. Does sparseData contain more?
        forAll(fullData, meshEdgeI)
        {
            const edge& e = mesh.edges()[meshEdgeI];

            const auto iter = sparseData.cfind(e);

            if (iter.good())
            {
                const point& sparsePt = iter.val();
                const point& fullPt = fullData[meshEdgeI];

                if (fullPt != sparsePt)
                {
                    FatalErrorInFunction
                        << "Extra edge:" << meshEdgeI
                        << " points:" << mesh.edges()[meshEdgeI]
                        << " full:" << fullPt
                        << " sparse:" << sparsePt
                        << exit(FatalError);
                }
            }
        }
    }
}


void testPointSync(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing point-wise data synchronisation." << endl;

    // Test position.

    {
        pointField syncedPoints(mesh.points());
        syncTools::syncPointPositions
        (
            mesh,
            syncedPoints,
            minMagSqrEqOp<point>(),
            point(GREAT, GREAT, GREAT)
        );

        forAll(syncedPoints, pointi)
        {
            if (mag(syncedPoints[pointi] - mesh.points()[pointi]) > SMALL)
            {
                FatalErrorInFunction
                    << "Point " << pointi
                    << " original location " << mesh.points()[pointi]
                    << " synced location " << syncedPoints[pointi]
                    << exit(FatalError);
            }
        }
    }

    // Test masterPoints

    {
        labelList nMasters(mesh.nPoints(), Zero);

        const bitSet isMasterPoint(syncTools::getMasterPoints(mesh));

        forAll(isMasterPoint, pointi)
        {
            if (isMasterPoint.test(pointi))
            {
                nMasters[pointi] = 1;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            nMasters,
            plusEqOp<label>(),
            0
        );

        forAll(nMasters, pointi)
        {
            if (nMasters[pointi] != 1)
            {
                WarningInFunction
                    << "Point " << pointi
                    << " original location " << mesh.points()[pointi]
                    << " has " << nMasters[pointi]
                    << " masters."
                    << endl;
            }
        }
    }
}


void testEdgeSync(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing edge-wise data synchronisation." << endl;

    const edgeList& edges = mesh.edges();

    // Test position.

    {
        pointField syncedMids(edges.size());
        forAll(syncedMids, edgeI)
        {
            syncedMids[edgeI] = edges[edgeI].centre(mesh.points());
        }
        syncTools::syncEdgePositions
        (
            mesh,
            syncedMids,
            minMagSqrEqOp<point>(),
            point(GREAT, GREAT, GREAT)
        );

        forAll(syncedMids, edgeI)
        {
            point eMid = edges[edgeI].centre(mesh.points());

            if (mag(syncedMids[edgeI] - eMid) > SMALL)
            {
                FatalErrorInFunction
                    << "Edge " << edgeI
                    << " original midpoint " << eMid
                    << " synced location " << syncedMids[edgeI]
                    << exit(FatalError);
            }
        }
    }

    // Test masterEdges

    {
        labelList nMasters(edges.size(), Zero);

        const bitSet isMasterEdge(syncTools::getMasterEdges(mesh));

        forAll(isMasterEdge, edgeI)
        {
            if (isMasterEdge.test(edgeI))
            {
                nMasters[edgeI] = 1;
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            nMasters,
            plusEqOp<label>(),
            0
        );

        forAll(nMasters, edgeI)
        {
            if (nMasters[edgeI] != 1)
            {
                const edge& e = edges[edgeI];
                WarningInFunction
                    << "Edge " << edgeI
                    << " at:" << mesh.points()[e[0]] << mesh.points()[e[1]]
                    << " has " << nMasters[edgeI]
                    << " masters."
                    << endl;
            }
        }
    }
}


typedef Pair<point> PointPair;
class edgePointCombineOp
{
public:
    void operator()(PointPair& x, const PointPair& y) const
    {
        if
        (
            (y[0] < x[0])
         || (y[0] == x[0] && y[1] < x[1])
        )
        {
            x = y;
        }
    }
};
class edgePointTransformOp
{
public:
    void operator()
    (
        const vectorTensorTransform& vt,
        const bool forward,
        List<PointPair>& fld
    ) const
    {
        pointField points0(fld.size());
        pointField points1(fld.size());
        forAll(fld, i)
        {
            points0[i] = fld[i].first();
            points1[i] = fld[i].second();
        }

        pointField points0New;
        pointField points1New;
        if (forward)
        {
            points0New = vt.transformPosition(points0);
            points1New = vt.transformPosition(points1);
        }
        else
        {
            points0New = vt.invTransformPosition(points0);
            points1New = vt.invTransformPosition(points1);
        }

        forAll(fld, i)
        {
            fld[i] = PointPair(points0New[i], points1New[i]);
        }
    }
};
class edgePointFlipOp
{
public:
    PointPair operator()(const PointPair& val) const
    {
        PointPair newVal(val);
        newVal.flip();
        return newVal;
    }
};
void testEdgeFlip2(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing edge-wise (oriented) data synchronisation." << endl;

    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    // Test position.

    List<PointPair> synEdgeEnds(edges.size());
    {
        forAll(synEdgeEnds, edgeI)
        {
            const edge& e = edges[edgeI];
            synEdgeEnds[edgeI] = PointPair
            (
                points[e[0]],
                points[e[1]]
            );
        }
    }

    // 1. Ignore flipping
    {
        List<PointPair> fld(synEdgeEnds);

        syncTools::syncEdgeList
        (
            mesh,
            fld,
            edgePointCombineOp(),
            PointPair(point::max, point::max),
            edgePointTransformOp(),
            identityOp()
        );

        forAll(fld, edgeI)
        {
            const edge& e = edges[edgeI];
            const PointPair edgeEnd
            (
                points[e[0]],
                points[e[1]]
            );

            const PointPair& sync = fld[edgeI];

            if
            (
                (mag(edgeEnd[0] - sync[0]) > SMALL)
             || (mag(edgeEnd[1] - sync[1]) > SMALL)
            )
            {
                WarningInFunction
                    << "Edge " << edgeI
                    << " original endpoints " << edgeEnd
                    << " synced endpoints " << sync
                    << endl;
            }
        }
    }

    // 2. Use flipping operator. Should produce no warnings
    {
        syncTools::syncEdgeList
        (
            mesh,
            synEdgeEnds,
            edgePointCombineOp(),
            PointPair(point::max, point::max),
            edgePointTransformOp(),
            edgePointFlipOp()
        );

        forAll(synEdgeEnds, edgeI)
        {
            const edge& e = edges[edgeI];
            const PointPair edgeEnd
            (
                points[e[0]],
                points[e[1]]
            );

            const PointPair& sync = synEdgeEnds[edgeI];

            if
            (
                (mag(edgeEnd[0] - sync[0]) > SMALL)
             || (mag(edgeEnd[1] - sync[1]) > SMALL)
            )
            {
                FatalErrorInFunction
                    << "Edge " << edgeI
                    << " original endpoints " << edgeEnd
                    << " synced endpoints " << sync
                    << exit(FatalError);
            }
        }
    }
}


void testEdgeFlip(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing edge-wise (oriented) data synchronisation."
        << endl;

    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    // Test vector.

    vectorField synEdgeVecs(edges.size());
    {
        forAll(synEdgeVecs, edgeI)
        {
            synEdgeVecs[edgeI] = edges[edgeI].unitVec(points);
        }
    }

    // Without flipping (should produce warnings)
    {
        vectorField fld(synEdgeVecs);
        // Ignore flipping
        syncTools::syncEdgeList
        (
            mesh,
            fld,
            minEqOp<vector>(),
            point::max
        );

        forAll(fld, edgeI)
        {
            const edge& e = edges[edgeI];
            const vector eVec(e.unitVec(points));

            if ((eVec & fld[edgeI]) < (1-SMALL))
            {
                WarningInFunction
                    << "Edge " << edgeI
                    << " at " << e.line(points)
                    << " original vector " << eVec
                    << " synced vector " << fld[edgeI]
                    << " diff:" << (eVec & fld[edgeI])
                    << endl;
            }
        }
    }


    // With consistent flipping. Should never produce difference
    {
        syncTools::syncEdgeList
        (
            mesh,
            synEdgeVecs,
            minMagSqrEqOp<vector>(),
            point::max,
            mapDistribute::transform(),
            flipOp()
        );

        forAll(synEdgeVecs, edgeI)
        {
            const edge& e = edges[edgeI];
            const vector eVec(e.unitVec(points));

            if ((eVec & synEdgeVecs[edgeI]) < (1-SMALL))
            {
                FatalErrorInFunction
                    << "Edge " << edgeI
                    << " at " << e.line(points)
                    << " original vector " << eVec
                    << " synced vector " << synEdgeVecs[edgeI]
                    << " diff:" << (eVec & synEdgeVecs[edgeI])
                    << exit(FatalError);
            }
        }
    }
}


class pointListOps
{
public:

    void operator()(pointList& lhs, const pointList& rhs) const
    {
        forAll(lhs, i)
        {
            point& l = lhs[i];
            const point& r = rhs[i];
            maxMagSqrEqOp<vector>()(l, r);
        }
    }

    void operator()
    (
        const vectorTensorTransform& vt,
        const bool forward,
        List<pointList>& fld
    ) const
    {
        if (forward)
        {
            for (auto& elems : fld)
            {
                for (auto& elem : elems)
                {
                    elem = vt.transformPosition(elem);
                }
            }
        }
        else
        {
            for (auto& elems : fld)
            {
                for (auto& elem : elems)
                {
                    elem = vt.invTransformPosition(elem);
                }
            }
        }
    }

    //- Transform patch-based field
    void operator()(const coupledPolyPatch& cpp, List<pointList>& fld) const
    {
        forAll(fld, facei)
        {
            pointList& pts = fld[facei];
            for (auto& pt : pts)
            {
                cpp.transformPosition(pt, facei);
            }
        }
    }
};


void testFaceSync(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing face-wise data synchronisation." << endl;

    // Test position.

    {
        pointField syncedFc(mesh.faceCentres());

        syncTools::syncFacePositions
        (
            mesh,
            syncedFc,
            maxMagSqrEqOp<point>()
        );

        forAll(syncedFc, facei)
        {
            if (mag(syncedFc[facei] - mesh.faceCentres()[facei]) > SMALL)
            {
                FatalErrorInFunction
                    << "Face " << facei
                    << " original centre " << mesh.faceCentres()[facei]
                    << " synced centre " << syncedFc[facei]
                    << exit(FatalError);
            }
        }
    }


    // Test non-contiguous data (uses streaming)

    {
        List<pointList> syncedFc(mesh.nFaces());

        const pointField& fcs = mesh.faceCentres();
        forAll(fcs, facei)
        {
            const point& fc = fcs[facei];
            syncedFc[facei].setSize(2, fc);
        }

        SubList<pointList> bndValues
        (
            syncedFc,
            mesh.nBoundaryFaces(),
            mesh.nInternalFaces()
        );
        syncTools::syncBoundaryFaceList
        (
            mesh,
            bndValues,
            pointListOps(), //does maxMagSqrEqOp<pointList>()
            pointListOps()  //transforms all
        );

        forAll(syncedFc, facei)
        {
            const point& fc = fcs[facei];
            const pointList& fld = syncedFc[facei];
            forAll(fld, i)
            {
                if (mag(fld[i] - fc) > SMALL)
                {
                    FatalErrorInFunction
                        << "Face " << facei
                        << " original centre " << fc
                        << " synced centre " << fld[i]
                        << exit(FatalError);
                }
            }
        }
    }


    // Test masterFaces

    {
        labelList nMasters(mesh.nFaces(), Zero);

        const bitSet isMasterFace(syncTools::getMasterFaces(mesh));

        forAll(isMasterFace, facei)
        {
            if (isMasterFace.test(facei))
            {
                nMasters[facei] = 1;
            }
        }

        syncTools::syncFaceList
        (
            mesh,
            nMasters,
            plusEqOp<label>()
        );

        forAll(nMasters, facei)
        {
            if (nMasters[facei] != 1)
            {
                FatalErrorInFunction
                    << "Face " << facei
                    << " centre " << mesh.faceCentres()[facei]
                    << " has " << nMasters[facei]
                    << " masters."
                    << exit(FatalError);
            }
        }
    }
}


// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"


    Random rndGen(5341*(Pstream::myProcNo()+1));


    // Face sync
    testFaceSync(mesh, rndGen);

    // Edge sync
    testEdgeSync(mesh, rndGen);

    // Edge sync and flip
    testEdgeFlip(mesh, rndGen);

    // Edge sync and flip of more complex structure
    testEdgeFlip2(mesh, rndGen);

    // Point sync
    testPointSync(mesh, rndGen);

    // PackedList synchronisation
    testPackedList(mesh, rndGen);

    // Sparse synchronisation
    testSparseData(mesh, rndGen);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
