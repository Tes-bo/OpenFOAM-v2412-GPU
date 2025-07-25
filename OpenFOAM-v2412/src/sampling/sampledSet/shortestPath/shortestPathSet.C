/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "shortestPathSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "topoDistanceData.H"
#include "addToRunTimeSelectionTable.H"
#include "FaceCellWave.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "volFields.H"
#include "OBJstream.H"
#include "PatchTools.H"
#include "foamVtkSurfaceWriter.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shortestPathSet, 0);
    addToRunTimeSelectionTable(sampledSet, shortestPathSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::shortestPathSet::findMinFace
(
    const polyMesh& mesh,
    const label cellI,
    const List<topoDistanceData<label>>& allFaceInfo,
    const bitSet& isLeakPoint,
    const bool distanceMode,
    const point& origin
)
{
    const cell& cFaces2 = mesh.cells()[cellI];

    // 1. Get topologically nearest face

    label minDist = labelMax;
    label minFaceI = -1;
    label nMin = 0;
    forAll(cFaces2, i)
    {
        label faceI = cFaces2[i];
        const topoDistanceData<label>& info = allFaceInfo[faceI];
        if (info.distance() < minDist)
        {
            minDist = info.distance();
            minFaceI = faceI;
            nMin = 1;
        }
        else if (info.distance() == minDist)
        {
            nMin++;
        }
    }

    if (nMin > 1)
    {
        // 2. Check all faces with minDist for minimum (or maximum)
        //    distance to origin
        if (distanceMode)
        {
            scalar minDist2 = ROOTVGREAT;
            forAll(cFaces2, i)
            {
                label faceI = cFaces2[i];
                if (allFaceInfo[faceI].distance() == minDist)
                {
                    scalar d2 = magSqr(mesh.faceCentres()[faceI]-origin);
                    if (d2 < minDist2)
                    {
                        minDist2 = d2;
                        minFaceI  = faceI;
                    }
                }
            }
        }
        else
        {
            // Avoid leak points i.e. preferentially stay away from the wall
            label minLeakPoints = labelMax;
            forAll(cFaces2, i)
            {
                label faceI = cFaces2[i];
                if (allFaceInfo[faceI].distance() == minDist)
                {
                    // Count number of leak points
                    label nLeak = 0;
                    {
                        const face& f = mesh.faces()[faceI];
                        forAll(f, fp)
                        {
                            if (isLeakPoint[f[fp]])
                            {
                                nLeak++;
                            }
                        }
                    }

                    if (nLeak < minLeakPoints)
                    {
                        minLeakPoints = nLeak;
                        minFaceI  = faceI;
                    }
                }
            }
        }
    }

    return minFaceI;
}


void Foam::shortestPathSet::calculateDistance
(
    const label iter,
    const polyMesh& mesh,
    const label cellI,

    List<topoDistanceData<label>>& allFaceInfo,
    List<topoDistanceData<label>>& allCellInfo
) const
{
    int dummyTrackData = 0;

    // Seed faces on cell1
    DynamicList<topoDistanceData<label>> faceDist;
    DynamicList<label> cFaces1;

    if (cellI != -1)
    {
        const labelList& cFaces = mesh.cells()[cellI];
        faceDist.reserve(cFaces.size());
        cFaces1.reserve(cFaces.size());

        for (label facei : cFaces)
        {
            if (!allFaceInfo[facei].valid(dummyTrackData))
            {
                cFaces1.append(facei);
                faceDist.append(topoDistanceData<label>(0, 123));
            }
        }
    }



    // Walk through face-cell wave till all cells are reached
    FaceCellWave
    <
        topoDistanceData<label>
    > wallDistCalc
    (
        mesh,
        cFaces1,
        faceDist,
        allFaceInfo,
        allCellInfo,
        mesh.globalData().nTotalCells()+1   // max iterations
    );


    if (debug)
    {
        const fvMesh& fm = refCast<const fvMesh>(mesh);

        //const_cast<fvMesh&>(fm).setInstance(fm.time().timeName());
        //fm.write();
        volScalarField fld
        (
            IOobject
            (
                "allCellInfo" + Foam::name(iter),
                fm.time().timeName(),
                fm,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fm,
            dimensionedScalar(dimless, Zero)
        );
        forAll(allCellInfo, celli)
        {
            fld[celli] = allCellInfo[celli].distance();
        }
        forAll(fld.boundaryField(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];
            SubList<topoDistanceData<label>> p(pp.patchSlice(allFaceInfo));
            scalarField pfld(fld.boundaryField()[patchi].size());
            forAll(pfld, i)
            {
                pfld[i] = 1.0*p[i].distance();
            }
            fld.boundaryFieldRef()[patchi] == pfld;
        }
        //Note: do not swap cell values so do not do
        //fld.correctBoundaryConditions();
        Pout<< "Writing distance field for initial cell " << cellI
            << " to " << fld.objectPath() << endl;
        fld.write();
    }
}


void Foam::shortestPathSet::sync
(
    const polyMesh& mesh,
    bitSet& isLeakFace,
    bitSet& isLeakPoint,
    const label celli,
    point& origin,
    bool& findMinDistance
) const
{
    syncTools::syncPointList
    (
        mesh,
        isLeakPoint,
        orEqOp<unsigned int>(),
        0u
    );
    syncTools::syncFaceList
    (
        mesh,
        isLeakFace,
        orEqOp<unsigned int>()
    );
    // Sync origin, findMinDistance
    {
        typedef Tuple2<label, Tuple2<point, bool>> ProcData;

        ProcData searchData
        (
            celli,
            Tuple2<point, bool>(origin, findMinDistance)
        );
        Pstream::combineReduce
        (
            searchData,
            [](ProcData& x, const ProcData& y)
            {
                if (y.first() != -1)
                {
                    x = y;
                }
            }
        );
        origin = searchData.second().first();
        findMinDistance = searchData.second().second();
    }
}


bool Foam::shortestPathSet::touchesWall
(
    const polyMesh& mesh,
    const label facei,

    bitSet& isLeakFace,
    const bitSet& isLeakPoint
) const
{
    // Check if facei touches leakPoint
    const face& f = mesh.faces()[facei];
    forAll(f, fp)
    {
        const label nextFp = f.fcIndex(fp);

        if (isLeakPoint[f[fp]] && isLeakPoint[f[nextFp]])   // edge on boundary
        //if (isLeakPoint[f[fp]])                           // point on boundary
        {
            if (isLeakFace.set(facei))
            {
                return true;
            }
        }
    }

    return false;
}


Foam::bitSet Foam::shortestPathSet::pathFaces
(
    const polyMesh& mesh,
    const bitSet& isLeakCell
) const
{
    // Calculates the list of faces inbetween leak cells.
    // Does not account for the fact that the cells might be from different
    // paths...

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& own = mesh.faceOwner();
    const labelList& nbr = mesh.faceNeighbour();

    // Get remote value of bitCell. Note: keep uncoupled boundary faces false.
    boolList nbrLeakCell(mesh.nBoundaryFaces(), false);
    {
        for (const polyPatch& pp : patches)
        {
            if (pp.coupled())
            {
                label bFacei = pp.start()-mesh.nInternalFaces();

                const labelUList& faceCells = pp.faceCells();

                for (const label celli : faceCells)
                {
                    nbrLeakCell[bFacei] = isLeakCell[celli];
                    ++bFacei;
                }
            }
        }

        syncTools::swapBoundaryFaceList
        (
            mesh,
            nbrLeakCell
        );
    }


    bitSet isLeakFace(mesh.nFaces());

    // Internal faces
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        if (isLeakCell[own[facei]] && isLeakCell[nbr[facei]])
        {
            isLeakFace.set(facei);
        }
    }
    // Boundary faces
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        if (isLeakCell[own[facei]] && nbrLeakCell[facei-mesh.nInternalFaces()])
        {
            isLeakFace.set(facei);
        }
    }
    return isLeakFace;
}


bool Foam::shortestPathSet::genSingleLeakPath
(
    const bool markLeakPath,
    const label iter,
    const polyMesh& mesh,
    const bitSet& isBlockedFace,
    const point& insidePoint,
    const label insideCelli,
    const point& outsidePoint,
    const label outsideCelli,

    // Generated sampling points
    const scalar trackLength,
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist,

    // State of current leak paths
    bitSet& isLeakCell,
    bitSet& isLeakFace,
    bitSet& isLeakPoint,

    // Work storage
    List<topoDistanceData<label>>& allFaceInfo,
    List<topoDistanceData<label>>& allCellInfo
) const
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const topoDistanceData<label> maxData(labelMax, labelMax);


    allFaceInfo.setSize(mesh.nFaces());
    allFaceInfo = topoDistanceData<label>();
    allCellInfo.setSize(mesh.nCells());
    allCellInfo = topoDistanceData<label>();

    // Mark blocked faces with high distance
    forAll(isBlockedFace, facei)
    {
        if (isBlockedFace[facei])
        {
            allFaceInfo[facei] = maxData;
        }
    }

    if (markLeakPath)
    {
        // Mark any previously found leak path. This marks all
        // faces of all cells on the path. This will make sure that
        // ultimately the path will go through another gap.
        forAll(isLeakCell, celli)
        {
            if (celli != insideCelli && celli != outsideCelli)
            {
                if (isLeakCell[celli])
                {
                    allCellInfo[celli] = maxData;
                    //- Mark all faces of the cell
                    //const cell& cFaces = mesh.cells()[celli];
                    //for (auto facei : cFaces)
                    //{
                    //    allFaceInfo[facei] = maxData;
                    //}
                }
            }
        }

        //- Mark only faces inbetween two leak cells
        bitSet isLeakCellWithout(isLeakCell);
        if (insideCelli != -1)
        {
            isLeakCellWithout.unset(insideCelli);
        }
        if (outsideCelli != -1)
        {
            isLeakCellWithout.unset(outsideCelli);
        }
        const bitSet isPathFace(pathFaces(mesh, isLeakCellWithout));
        forAll(isPathFace, facei)
        {
            if (isPathFace[facei])
            {
                allFaceInfo[facei] = maxData;
            }
        }
    }

    // Mark any previously found leak faces. These are faces that
    // are (point- or edge-)connected to an existing boundary.
    forAll(isLeakFace, facei)
    {
        if (isLeakFace[facei])
        {
            allFaceInfo[facei] = maxData;
        }
    }


    // Pass1: Get distance to insideCelli

    calculateDistance(iter, mesh, insideCelli, allFaceInfo, allCellInfo);



    // Pass2: walk from outside points backwards. Note: could be done
    //        using FaceCellWave as well but is overly complex since
    //        does not allow logic comparing all faces of a cell.

    bool targetFound = false;
    if (outsideCelli != -1)
    {
        int dummyTrackData;
        targetFound = allCellInfo[outsideCelli].valid(dummyTrackData);
        if (iter == 0 && !targetFound)
        {
            WarningInFunction
                << "Point " << outsidePoint
                << " not reachable by walk from " << insidePoint
                << ". Probably mesh has island/regions."
                << " Skipped route detection." << endl;
        }
    }

    if (!returnReduceOr(targetFound))
    {
        //Pout<< "now :"
        //    << " nLeakCell:"
        //    << returnReduce(isLeakCell.count(), sumOp<label>())
        //    << " nLeakFace:"
        //    << returnReduce(isLeakFace.count(), sumOp<label>())
        //    << " nLeakPoint:"
        //    << returnReduce(isLeakPoint.count(), sumOp<label>())
        //    << endl;

        return false;
    }


    // Start with given target cell and walk back
    // If point happens to be on multiple processors, random pick
    label frontCellI = outsideCelli;
    point origin(outsidePoint);
    bool findMinDistance = true;

    while (true)
    {
        // We have a single face front (frontFaceI). Walk from face to cell
        // to face etc until we reach the destination cell.

        label frontFaceI = -1;

        // Work within same processor
        if (frontCellI != -1)
        {
            // Find face with lowest distance from seed. In below figure the
            // seed cell is marked with distance 0. It is surrounded by faces
            // and cells with distance 1. The layer outside is marked with
            // distance 2 etc etc.
            //
            //   x  |  x  2  1  2  2  |  x  |  x
            //  --- + --- + -1- + -2- + --- + ---
            //   x  |  1  1  0  1  1  |  x  |  x
            //  --- + --- + -1- + -2- + --- + ---
            //   x  |  x  2  1  2  2  3  3  4  4
            //  --- + --- + --- + -3- + -4- + -5-
            //   x  |  x  3  2  3  3  4  4  5  5
            //
            // e.g. if we start backwards search from cell with value = 4,
            // we have neighbour faces 4, 4, 5, 5. Choose 4 (least distance
            // to seed) and continue...

            frontFaceI = findMinFace
            (
                mesh,
                frontCellI,
                allFaceInfo,
                isLeakPoint,
                findMinDistance,    // mode : find min or find max
                origin
            );


            // Loop until we hit a boundary face
            bitSet isNewLeakPoint(isLeakPoint);
            while (mesh.isInternalFace(frontFaceI))
            {
                if (isBlockedFace.size() && isBlockedFace[frontFaceI])
                {
                    // This should not happen since we never walk through
                    // a blocked face. However ...
                    // Pout<< " Found blocked face" << endl;
                    frontCellI = -1;
                    break;
                }

                // Step to neighbouring cell
                label nbrCellI = mesh.faceOwner()[frontFaceI];
                if (nbrCellI == frontCellI)
                {
                    nbrCellI = mesh.faceNeighbour()[frontFaceI];
                }

                if (nbrCellI == insideCelli)
                {
                    // Reached destination. This is the normal exit.
                    frontCellI = -1;
                    break;
                }

                frontCellI = nbrCellI;

                // Pick best face on cell
                frontFaceI = findMinFace
                (
                    mesh,
                    frontCellI,
                    allFaceInfo,
                    isLeakPoint,
                    findMinDistance,
                    origin
                );

                const topoDistanceData<label>& cInfo = allCellInfo[frontCellI];
                const topoDistanceData<label>& fInfo = allFaceInfo[frontFaceI];

                if (fInfo.distance() <= cInfo.distance())
                {
                    // Found valid next cell,face. Mark on path
                    samplingPts.append(mesh.cellCentres()[frontCellI]);
                    samplingCells.append(frontCellI);
                    samplingFaces.append(-1);
                    samplingSegments.append(iter);
                    samplingCurveDist.append
                    (
                        trackLength+cInfo.distance()
                    );
                    isLeakCell.set(frontCellI);

                    // Check if front face uses a boundary point.
                    if
                    (
                        touchesWall
                        (
                            mesh,
                            frontFaceI,
                            isLeakFace,
                            isLeakPoint
                        )
                    )
                    {
                        //Pout<< "** added leak face:" << frontFaceI << " at:"
                        //<< mesh.faceCentres()[frontFaceI] << endl;
                        isNewLeakPoint.set(mesh.faces()[frontFaceI]);
                        origin = mesh.faceCentres()[frontFaceI];
                        findMinDistance = false;
                    }
                }
            }
            isLeakPoint.transfer(isNewLeakPoint);
        }

        // Situation 1: we found the destination cell (do nothing),
        // frontCellI is -1 on all processors
        if (returnReduceAnd(frontCellI == -1))
        {
            //Pout<< "now :"
            //    << " nLeakCell:"
            //    << returnReduce(isLeakCell.count(), sumOp<label>())
            //    << " nLeakFace:"
            //    << returnReduce(isLeakFace.count(), sumOp<label>())
            //    << " nLeakPoint:"
            //    << returnReduce(isLeakPoint.count(), sumOp<label>())
            //    << endl;

            break;
        }

        // Situation 2: we're on a coupled patch and might need to
        //              switch processor/cell. We need to transfer:
        //              -frontface -frontdistance -leak points/faces
        boolList isFront(mesh.nBoundaryFaces(), false);

        if (frontFaceI != -1)
        {
            isFront[frontFaceI-mesh.nInternalFaces()] = true;
        }
        syncTools::swapBoundaryFaceList(mesh, isFront);

        label minCellDistance = labelMax;
        if (frontCellI != -1)
        {
            minCellDistance = allCellInfo[frontCellI].distance();
        }
        reduce(minCellDistance, minOp<label>());

        // Sync all leak data
        sync
        (
            mesh,
            isLeakFace,
            isLeakPoint,
            frontCellI,
            origin,
            findMinDistance
        );


        // Bit tricky:
        // - processor might get frontFaceI/frontCellI in through multiple faces
        //   (even through different patches?)
        // - so loop over all (coupled) patches and pick the best frontCellI

        const label oldFrontFaceI = frontFaceI;
        frontCellI = -1;
        frontFaceI = -1;
        forAll(pbm, patchI)
        {
            const polyPatch& pp = pbm[patchI];
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                if
                (
                    oldFrontFaceI == -1
                 && isFront[faceI-mesh.nInternalFaces()]
                 && (isBlockedFace.empty() || !isBlockedFace[faceI])
                )
                {
                    frontFaceI = faceI;
                    frontCellI = pp.faceCells()[i];
                    break;
                }
            }

            if
            (
                frontCellI != -1
             && allCellInfo[frontCellI].distance() < minCellDistance
            )
            {
                const topoDistanceData<label>& cInfo = allCellInfo[frontCellI];

                samplingPts.append(mesh.cellCentres()[frontCellI]);
                samplingCells.append(frontCellI);
                samplingFaces.append(-1);
                samplingSegments.append(iter);
                samplingCurveDist.append
                (
                    trackLength+cInfo.distance()
                );
                isLeakCell.set(frontCellI);

                // Check if frontFacei touches leakPoint
                if
                (
                    touchesWall
                    (
                        mesh,
                        frontFaceI,
                        isLeakFace,
                        isLeakPoint
                    )
                )
                {
                    //Pout<< "** added leak BOUNDARY face:" << frontFaceI
                    // << " at:" << mesh.faceCentres()[frontFaceI] << endl;
                    isLeakPoint.set(mesh.faces()[frontFaceI]);
                    origin = mesh.faceCentres()[frontFaceI];
                    findMinDistance = false;
                }

                break;
            }

            // When seed cell is isolated by processor boundaries
            if (insideCelli != -1 && frontCellI == insideCelli)
            {
                // Pout<< " Found connection seed cell!" << endl;
                frontCellI = -1;
                break;
            }
        }

        // Sync all leak data
        sync
        (
            mesh,
            isLeakFace,
            isLeakPoint,
            frontCellI,
            origin,
            findMinDistance
        );
    }

    return true;
}


// Clean up leak faces (erode open edges). These are leak faces which are
// not connected to another leak face or leak point. Parallel consistent.
// Returns overall number of faces deselected.
Foam::label Foam::shortestPathSet::erodeFaceSet
(
    const polyMesh& mesh,
    const bitSet& isBlockedPoint,
    bitSet& isLeakFace
) const
{
    if
    (
        (isBlockedPoint.size() != mesh.nPoints())
     || (isLeakFace.size() != mesh.nFaces())
    )
    {
        FatalErrorInFunction << "Problem :"
            << " isBlockedPoint:" << isBlockedPoint.size()
            << " isLeakFace:" << isLeakFace.size()
            << exit(FatalError);
    }

    const globalMeshData& globalData = mesh.globalData();
    const mapDistribute& map = globalData.globalEdgeSlavesMap();
    const indirectPrimitivePatch& cpp = mesh.globalData().coupledPatch();


    label nTotalEroded = 0;

    while (true)
    {
        bitSet newIsLeakFace(isLeakFace);

        // Get number of edges

        const labelList meshFaceIDs(isLeakFace.toc());
        const uindirectPrimitivePatch pp
        (
            UIndirectList<face>(mesh.faces(), meshFaceIDs),
            mesh.points()
        );

        // Count number of faces per edge
        const labelListList& edgeFaces = pp.edgeFaces();
        labelList nEdgeFaces(edgeFaces.size());
        forAll(edgeFaces, edgei)
        {
            nEdgeFaces[edgei] = edgeFaces[edgei].size();
        }

        // Match pp edges to coupled edges
        labelList patchEdges;
        labelList coupledEdges;
        bitSet sameEdgeOrientation;
        PatchTools::matchEdges
        (
            pp,
            cpp,
            patchEdges,
            coupledEdges,
            sameEdgeOrientation
        );


        // Convert patch-edge data into cpp-edge data
        labelList coupledNEdgeFaces(map.constructSize(), Zero);
        UIndirectList<label>(coupledNEdgeFaces, coupledEdges) =
            UIndirectList<label>(nEdgeFaces, patchEdges);

        // Synchronise
        globalData.syncData
        (
            coupledNEdgeFaces,
            globalData.globalEdgeSlaves(),
            globalData.globalEdgeTransformedSlaves(),
            map,
            plusEqOp<label>()
        );

        // Convert back from cpp-edge to patch-edge
        UIndirectList<label>(nEdgeFaces, patchEdges) =
            UIndirectList<label>(coupledNEdgeFaces, coupledEdges);

        // Remove any open edges
        label nEroded = 0;
        forAll(nEdgeFaces, edgei)
        {
            if (nEdgeFaces[edgei] == 1)
            {
                const edge& e = pp.edges()[edgei];
                const label mp0 = pp.meshPoints()[e[0]];
                const label mp1 = pp.meshPoints()[e[1]];

                if (!isBlockedPoint[mp0] || !isBlockedPoint[mp1])
                {
                    // Edge is not on wall so is open
                    const label patchFacei = edgeFaces[edgei][0];
                    const label meshFacei = pp.addressing()[patchFacei];
                    //Pout<< "Eroding face:" << meshFacei
                    //    << " at:" << mesh.faceCentres()[meshFacei]
                    //    << " since has open edge:" << mesh.points()[mp0]
                    //    << mesh.points()[mp1] << endl;

                    if (newIsLeakFace.unset(meshFacei))
                    {
                        nEroded++;
                    }
                }
            }
        }

        reduce(nEroded, sumOp<label>());
        nTotalEroded += nEroded;

        if (nEroded == 0)
        {
            break;
        }

        // Make sure we make the same decision on both sides
        syncTools::syncFaceList
        (
            mesh,
            newIsLeakFace,
            orEqOp<unsigned int>()
        );

        isLeakFace = std::move(newIsLeakFace);
    }

    return nTotalEroded;
}


void Foam::shortestPathSet::genSamples
(
    const bool addLeakPath,
    const label maxIter,
    const polyMesh& mesh,
    const bitSet& isBoundaryFace,
    const point& insidePoint,
    const label insideCelli,
    const point& outsidePoint,

    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist,
    bitSet& isLeakCell,
    bitSet& isLeakFace,
    bitSet& isLeakPoint
) const
{
    // Mark all paths needed to close a single combination of insidePoint,
    // outsidePoint. The output is
    // - isLeakCell : is cell along a leak path
    // - isLeakFace : is face along a leak path (so inbetween two leakCells)
    // - isLeakPoint : is point on a leakFace


    const topoDistanceData<label> maxData(labelMax, labelMax);

    // Get the target point
    const label outsideCelli = mesh.findCell(outsidePoint);

    // Maintain overall track length. Used to make curve distance continuous.
    scalar trackLength = 0;

    List<topoDistanceData<label>> allFaceInfo(mesh.nFaces());
    List<topoDistanceData<label>> allCellInfo(mesh.nCells());


    // Boundary face + additional temporary blocks (to force leakpath to
    // outside)
    autoPtr<bitSet> isBlockedFace;

    label iter;
    bool markLeakPath = false;


    for (iter = 0; iter < maxIter; iter++)
    {
        const label nOldLeakFaces = returnReduce
        (
            isLeakFace.count(),
            sumOp<label>()
        );
        const label oldNSamplingPts(samplingPts.size());

        bool foundPath = genSingleLeakPath
        (
            markLeakPath,
            iter,
            mesh,
            (isBlockedFace ? *isBlockedFace : isBoundaryFace),
            insidePoint,
            insideCelli,
            outsidePoint,
            outsideCelli,

            // Generated sampling points
            trackLength,
            samplingPts,
            samplingCells,
            samplingFaces,
            samplingSegments,
            samplingCurveDist,

            // State of current leak paths
            isLeakCell,
            isLeakFace,
            isLeakPoint,

            // Work storage
            allFaceInfo,
            allCellInfo
        );

        // Recalculate the overall trackLength
        trackLength = returnReduce
        (
            (
                samplingCurveDist.size()
              ? samplingCurveDist.last()
              : 0
            ),
            maxOp<scalar>()
        );

        const label nLeakFaces = returnReduce
        (
            isLeakFace.count(),
            sumOp<label>()
        );

        if (!foundPath && !markLeakPath)
        {
            // In mark-boundary-face-only mode and already fully blocked the
            // path to outsideCell so we're good
            break;
        }


        if (nLeakFaces > nOldLeakFaces)
        {
            // Normal operation: walking has closed some wall-connected faces
            // If previous iteration was markLeakPath-mode make sure to revert
            // to normal operation (i.e. face marked in isLeakFace)
            isBlockedFace.reset(nullptr);
            markLeakPath = false;
        }
        else
        {
            // Did not mark any additional faces/points. Save current state
            // and add faces/points on leakpath to force next walk
            // to pass outside of leakpath. This is done until the leakpath
            // 'touchesWall' (i.e. edge connected to an original boundary face)

            if (markLeakPath && !foundPath)
            {
                // Is marking all faces on all paths and no additional path
                // found. Also nothing new marked (in isLeakFace) since
                // nLeakFaces == nOldLeakFaces) so we're
                // giving up. Convert all path faces into leak faces
                //Pout<< "** giving up" << endl;
                break;
            }


            // Revert to boundaryFaces only
            if (!isBlockedFace)
            {
                //Pout<< "** Starting from original boundary faces." << endl;
                isBlockedFace.reset(new bitSet(isBoundaryFace));
            }

            markLeakPath = true;


            if (debug & 2)
            {
                const pointField leakCcs(mesh.cellCentres(), isLeakCell.toc());
                mkDir(mesh.time().timePath());
                OBJstream str
                (
                    mesh.time().timePath()
                   /"isLeakCell" + Foam::name(iter) + ".obj"
                );
                Pout<< "Writing new isLeakCell to " << str.name() << endl;
                str.write(leakCcs);
            }
            if (debug & 2)
            {
                OBJstream str
                (
                    mesh.time().timePath()
                   /"leakPath" + Foam::name(iter) + ".obj"
                );
                Pout<< "Writing new leak-path to " << str.name() << endl;
                for
                (
                    label samplei = oldNSamplingPts+1;
                    samplei < samplingPts.size();
                    samplei++
                )
                {
                    Pout<< "    passing through cell "
                        << samplingCells[samplei]
                        << " at:" << mesh.cellCentres()[samplingCells[samplei]]
                        << " distance:" << samplingCurveDist[samplei]
                        << endl;

                    str.writeLine
                    (
                        samplingPts[samplei-1],
                        samplingPts[samplei]
                    );
                }
            }
        }
    }

    if (debug)
    {
        const fvMesh& fm = refCast<const fvMesh>(mesh);

        const_cast<fvMesh&>(fm).setInstance(fm.time().timeName());
        volScalarField fld
        (
            IOobject
            (
                "isLeakCell",
                fm.time().timeName(),
                fm,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fm,
            dimensionedScalar(dimless, Zero)
        );
        for (const label celli : isLeakCell)
        {
            fld[celli] = scalar(1);
        }
        fld.correctBoundaryConditions();
        fld.write();
    }

    if (maxIter > 1 && iter == maxIter)
    {
        WarningInFunction << "Did not manage to close gap using " << iter
            << " leak paths" << nl << "This can cause problems when using the"
            << " paths to close leaks" << endl;
    }
}


void Foam::shortestPathSet::genSamples
(
    const bool markLeakPath,
    const label maxIter,
    const polyMesh& mesh,
    const labelUList& wallPatches,
    const bitSet& isBlockedFace
)
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    // Seed faces and points on 'real' boundary
    bitSet isBlockedPoint(mesh.nPoints());
    {
        // Real boundaries
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        for (label patchi : wallPatches)
        {
            const polyPatch& pp = pbm[patchi];
            forAll(pp, i)
            {
                isBlockedPoint.set(pp[i]);
            }
        }

        // Meshed boundaries
        forAll(isBlockedFace, facei)
        {
            if (isBlockedFace[facei])
            {
                isBlockedPoint.set(mesh.faces()[facei]);
            }
        }

        syncTools::syncPointList
        (
            mesh,
            isBlockedPoint,
            orEqOp<unsigned int>(),
            0u
        );

        if (debug)
        {
            mkDir(mesh.time().timePath());
            OBJstream str(mesh.time().timePath()/"isBlockedPoint.obj");
            for (const auto& pointi : isBlockedPoint)
            {
                str.write(mesh.points()[pointi]);
            }
            Pout<< "Writing " << str.nVertices()
                << " points to " << str.name() << endl;
        }
    }


    bitSet isLeakPoint(isBlockedPoint);
    // Newly closed faces
    bitSet isLeakFace(mesh.nFaces());
    // All cells along leak paths
    bitSet isLeakCell(mesh.nCells());

    label prevSegmenti = 0;
    scalar prevDistance = 0.0;

    for (auto insidePoint : insidePoints_)
    {
        const label insideCelli = mesh.findCell(insidePoint);

        for (auto outsidePoint : outsidePoints_)
        {
            const label nOldSamples = samplingSegments.size();

            // Append any leak path to sampling*
            genSamples
            (
                markLeakPath,
                maxIter,
                mesh,
                isBlockedFace,
                insidePoint,
                insideCelli,
                outsidePoint,

                samplingPts,
                samplingCells,
                samplingFaces,
                samplingSegments,
                samplingCurveDist,
                isLeakCell,
                isLeakFace,
                isLeakPoint
            );

            // Make segment, distance consecutive
            label maxSegment = 0;
            scalar maxDistance = 0.0;
            for (label i = nOldSamples; i < samplingSegments.size(); ++i)
            {
                samplingSegments[i] += prevSegmenti;
                maxSegment = max(maxSegment, samplingSegments[i]);
                samplingCurveDist[i] += prevDistance;
                maxDistance = max(maxDistance, samplingCurveDist[i]);
            }
            prevSegmenti = returnReduce(maxSegment, maxOp<label>());
            prevDistance = returnReduce(maxDistance, maxOp<scalar>());
        }
    }

    // Clean up leak faces (erode open edges). These are leak faces which are
    // not connected to another leak face or leak point
    erodeFaceSet(mesh, isBlockedPoint, isLeakFace);

    leakFaces_ = isLeakFace.sortedToc();


    if (debug)
    {
        const faceList leakFaces(mesh.faces(), leakFaces_);

        mkDir(mesh.time().timePath());
        OBJstream str(mesh.time().timePath()/"isLeakFace.obj");
        str.write(leakFaces, mesh.points(), false);
        Pout<< "Writing " << leakFaces.size()
            << " faces to " << str.name() << endl;
    }


    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    // Move into *this
    setSamples
    (
        std::move(samplingPts),
        std::move(samplingCells),
        std::move(samplingFaces),
        std::move(samplingSegments),
        std::move(samplingCurveDist)
    );

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shortestPathSet::shortestPathSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const bool markLeakPath,
    const label maxIter,
    const labelUList& wallPatches,
    const pointField& insidePoints,
    const pointField& outsidePoints,
    const boolList& isBlockedFace
)
:
    sampledSet(name, mesh, searchEngine, axis),
    insidePoints_(insidePoints),
    outsidePoints_(outsidePoints)
{
    if (debug)
    {
        fileName outputDir
        (
            mesh.time().globalPath()
          / functionObject::outputPrefix
          / mesh.pointsInstance()
        );
        outputDir.clean();  // Remove unneeded ".."

        Info<< "shortestPathSet : Writing blocked faces to "
            << outputDir << endl;

        const indirectPrimitivePatch setPatch
        (
            IndirectList<face>
            (
                mesh.faces(),
                findIndices(isBlockedFace, true)
            ),
            mesh.points()
        );

        if (Pstream::parRun())
        {
            // Topological merge
            labelList pointToGlobal;
            labelList uniqueMeshPointLabels;
            autoPtr<globalIndex> globalPoints;
            autoPtr<globalIndex> globalFaces;
            faceList mergedFaces;
            pointField mergedPoints;
            Foam::PatchTools::gatherAndMerge
            (
                mesh,
                setPatch.localFaces(),
                setPatch.meshPoints(),
                setPatch.meshPointMap(),

                pointToGlobal,
                uniqueMeshPointLabels,
                globalPoints,
                globalFaces,

                mergedFaces,
                mergedPoints
            );

            // Write
            if (Pstream::master())
            {
                vtk::surfaceWriter writer
                (
                    mergedPoints,
                    mergedFaces,
                    (outputDir / "blockedFace"),
                    false  // serial only - already merged
                );

                writer.writeGeometry();
            }
        }
        else
        {
            vtk::surfaceWriter writer
            (
                setPatch.localPoints(),
                setPatch.localFaces(),
                (outputDir / "blockedFace"),
                false  // serial only - redundant
            );

            writer.writeGeometry();
        }
    }

    genSamples
    (
        markLeakPath,
        maxIter,
        mesh,
        wallPatches,
        bitSet(isBlockedFace)
    );
}


Foam::shortestPathSet::shortestPathSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    insidePoints_(dict.get<pointField>("insidePoints")),
    outsidePoints_(dict.get<pointField>("outsidePoints"))
{
    const label maxIter(dict.getOrDefault<label>("maxIter", 50));
    const bool markLeakPath(dict.getOrDefault("markLeakPath", true));

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    DynamicList<label> wallPatches(pbm.size());
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (!pp.coupled() && !isA<emptyPolyPatch>(pp))
        {
            wallPatches.append(patchi);
        }
    }

    genSamples(markLeakPath, maxIter, mesh, wallPatches, bitSet());
}


// ************************************************************************* //
