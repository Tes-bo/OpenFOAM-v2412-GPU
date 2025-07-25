/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "syncTools.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "globalMeshData.H"
#include "transformList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class CombineOp>
void Foam::syncTools::combine
(
    Map<T>& pointValues,
    const CombineOp& cop,
    const label index,
    const T& val
)
{
    auto iter = pointValues.find(index);

    if (iter.good())
    {
        cop(iter.val(), val);
    }
    else
    {
        pointValues.insert(index, val);
    }
}


template<class T, class CombineOp>
void Foam::syncTools::combine
(
    EdgeMap<T>& edgeValues,
    const CombineOp& cop,
    const edge& index,
    const T& val
)
{
    auto iter = edgeValues.find(index);

    if (iter.good())
    {
        cop(iter.val(), val);
    }
    else
    {
        edgeValues.insert(index, val);
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointMap
(
    const polyMesh& mesh,
    Map<T>& pointValues,        // from mesh point label to value
    const CombineOp& cop,
    const TransformOp& top
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Make sure we use unique message tag
    const int oldTag = UPstream::incrMsgType();

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    // Values on shared points. Keyed on global shared index.
    Map<T> sharedPointValues;

    if (pd.nGlobalPoints() > 0)
    {
        // meshPoint per local index
        const labelList& sharedPtLabels = pd.sharedPointLabels();

        // global shared index per local index
        const labelList& sharedPtAddr = pd.sharedPointAddr();

        sharedPointValues.reserve(sharedPtAddr.size());

        // Fill my entries in the shared points
        forAll(sharedPtLabels, i)
        {
            const auto fnd = pointValues.cfind(sharedPtLabels[i]);

            if (fnd.good())
            {
                combine
                (
                    sharedPointValues,
                    cop,
                    sharedPtAddr[i],    // index
                    fnd.val()           // value
                );
            }
        }
    }


    if (UPstream::parRun())
    {
        // Presize according to number of processor patches
        // (global topology information may not yet be available...)
        DynamicList<label> neighbProcs(patches.nProcessorPatches());
        PstreamBuffers pBufs;

        // Reduce communication by only sending non-zero data,
        // but with multiply-connected processor/processor
        // (eg, processorCyclic) also need to send zero information
        // to keep things synchronised

        // Initialize for registerSend() bookkeeping
        pBufs.initRegisterSend();


        // Sample and send.

        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.nPoints())
            {
                const auto& procPatch = *ppp;
                const label nbrProci = procPatch.neighbProcNo();

                // Get data per patchPoint in neighbouring point numbers.
                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                // Extract local values. Create map from nbrPoint to value.
                // Note: how small initial size?
                Map<T> patchInfo(meshPts.size() / 20);

                forAll(meshPts, i)
                {
                    const auto iter = pointValues.cfind(meshPts[i]);

                    if (iter.good())
                    {
                        patchInfo.insert(nbrPts[i], iter.val());
                    }
                }

                // Neighbour connectivity
                neighbProcs.push_uniq(nbrProci);

                // Send to neighbour
                {
                    UOPstream toNbr(nbrProci, pBufs);
                    toNbr << patchInfo;

                    // Record if send is required (data are non-zero)
                    pBufs.registerSend(nbrProci, !patchInfo.empty());
                }
            }
        }

        // Limit exchange to involved procs.
        // - automatically discards unnecessary (unregistered) sends
        pBufs.finishedNeighbourSends(neighbProcs);


        // Receive and combine.
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.nPoints())
            {
                const auto& procPatch = *ppp;
                const label nbrProci = procPatch.neighbProcNo();

                if (!pBufs.recvDataCount(nbrProci))
                {
                    // Nothing to receive
                    continue;
                }

                Map<T> nbrPatchInfo;
                {
                    UIPstream fromNbr(nbrProci, pBufs);
                    fromNbr >> nbrPatchInfo;
                }

                // Transform
                top(procPatch, nbrPatchInfo);

                const labelList& meshPts = procPatch.meshPoints();

                // Only update those values which come from neighbour
                forAllConstIters(nbrPatchInfo, nbrIter)
                {
                    combine
                    (
                        pointValues,
                        cop,
                        meshPts[nbrIter.key()],
                        nbrIter.val()
                    );
                }
            }
        }
    }

    // Do the cyclics.
    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            // Owner does all.

            const cyclicPolyPatch& cycPatch = *cpp;
            const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPtsA = cycPatch.meshPoints();
            const labelList& meshPtsB = nbrPatch.meshPoints();

            // Extract local values. Create map from coupled-edge to value.
            Map<T> half0Values(meshPtsA.size() / 20);
            Map<T> half1Values(half0Values.size());

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                const auto point0Fnd = pointValues.cfind(meshPtsA[e[0]]);

                if (point0Fnd.good())
                {
                    half0Values.insert(i, *point0Fnd);
                }

                const auto point1Fnd = pointValues.cfind(meshPtsB[e[1]]);

                if (point1Fnd.good())
                {
                    half1Values.insert(i, *point1Fnd);
                }
            }

            // Transform to receiving side
            top(cycPatch, half1Values);
            top(nbrPatch, half0Values);

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                const auto half0Fnd = half0Values.cfind(i);

                if (half0Fnd.good())
                {
                    combine
                    (
                        pointValues,
                        cop,
                        meshPtsB[e[1]],
                        *half0Fnd
                    );
                }

                const auto half1Fnd = half1Values.cfind(i);

                if (half1Fnd.good())
                {
                    combine
                    (
                        pointValues,
                        cop,
                        meshPtsA[e[0]],
                        *half1Fnd
                    );
                }
            }
        }
    }

    // Synchronize multiple shared points.
    if (pd.nGlobalPoints() > 0)
    {
        // meshPoint per local index
        const labelList& sharedPtLabels = pd.sharedPointLabels();
        // global shared index per local index
        const labelList& sharedPtAddr = pd.sharedPointAddr();

        // Reduce on master.

        if (UPstream::parRun())
        {
            if (UPstream::master())
            {
                // Receive the edges using shared points from other procs
                for (const int proci : UPstream::subProcs())
                {
                    Map<T> nbrValues;
                    IPstream::recv(nbrValues, proci);

                    // Merge neighbouring values with my values
                    forAllConstIters(nbrValues, iter)
                    {
                        combine
                        (
                            sharedPointValues,
                            cop,
                            iter.key(),     // edge
                            iter.val()      // value
                        );
                    }
                }
            }
            else
            {
                // Send to master
                OPstream::send(sharedPointValues, UPstream::masterNo());
            }

            // Broadcast: send merged values to all
            Pstream::broadcast(sharedPointValues);
        }

        // Merge sharedPointValues (keyed on sharedPointAddr) into
        // pointValues (keyed on mesh points).

        // Map from global shared index to meshpoint
        Map<label> sharedToMeshPoint(sharedPtAddr, sharedPtLabels);

        forAllConstIters(sharedToMeshPoint, iter)
        {
            // Do I have a value for my shared point
            const auto sharedFnd = sharedPointValues.cfind(iter.key());

            if (sharedFnd.good())
            {
                pointValues.set(iter.val(), sharedFnd.val());
            }
        }
    }

    // Reset tag
    UPstream::msgType(oldTag);
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncEdgeMap
(
    const polyMesh& mesh,
    EdgeMap<T>& edgeValues,
    const CombineOp& cop,
    const TransformOp& top
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Do synchronisation without constructing globalEdge addressing
    // (since this constructs mesh edge addressing)

    // Make sure we use unique message tag
    const int oldTag = UPstream::incrMsgType();

    // Swap proc patch info
    // ~~~~~~~~~~~~~~~~~~~~

    if (UPstream::parRun())
    {
        // Presize according to number of processor patches
        // (global topology information may not yet be available...)
        DynamicList<label> neighbProcs(patches.nProcessorPatches());
        PstreamBuffers pBufs;

        // Reduce communication by only sending non-zero data,
        // but with multiply-connected processor/processor
        // (eg, processorCyclic) also need to send zero information
        // to keep things synchronised

        // Initialize for registerSend() bookkeeping
        pBufs.initRegisterSend();


        // Sample and send.

        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.nEdges())
            {
                const auto& procPatch = *ppp;
                const label nbrProci = procPatch.neighbProcNo();

                // Get data per patch edge in neighbouring edge.
                const edgeList& edges = procPatch.edges();
                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                EdgeMap<T> patchInfo(edges.size() / 20);

                for (const edge& e : edges)
                {
                    const edge meshEdge(meshPts, e);

                    const auto iter = edgeValues.cfind(meshEdge);

                    if (iter.good())
                    {
                        const edge nbrEdge(nbrPts, e);
                        patchInfo.insert(nbrEdge, iter.val());
                    }
                }


                // Neighbour connectivity
                neighbProcs.push_uniq(nbrProci);

                // Send to neighbour
                {
                    UOPstream toNbr(nbrProci, pBufs);
                    toNbr << patchInfo;

                    // Record if send is required (data are non-zero)
                    pBufs.registerSend(nbrProci, !patchInfo.empty());
                }
            }
        }

        // Limit exchange to involved procs
        // - automatically discards unnecessary (unregistered) sends
        pBufs.finishedNeighbourSends(neighbProcs);


        // Receive and combine.
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.nEdges())
            {
                const auto& procPatch = *ppp;
                const label nbrProci = procPatch.neighbProcNo();

                if (!pBufs.recvDataCount(nbrProci))
                {
                    // Nothing to receive
                    continue;
                }

                EdgeMap<T> nbrPatchInfo;
                {
                    UIPstream fromNbr(nbrProci, pBufs);
                    fromNbr >> nbrPatchInfo;
                }

                // Apply transform to convert to this side properties.
                top(procPatch, nbrPatchInfo);


                // Only update those values which come from neighbour
                const labelList& meshPts = procPatch.meshPoints();

                forAllConstIters(nbrPatchInfo, nbrIter)
                {
                    const edge& e = nbrIter.key();
                    const edge meshEdge(meshPts, e);

                    combine
                    (
                        edgeValues,
                        cop,
                        meshEdge,           // edge
                        nbrIter.val()       // value
                    );
                }
            }
        }
    }


    // Swap cyclic info
    // ~~~~~~~~~~~~~~~~

    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            // Owner does all.

            const cyclicPolyPatch& cycPatch = *cpp;
            const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();

            const edgeList& coupledEdges = cycPatch.coupledEdges();

            const labelList& meshPtsA = cycPatch.meshPoints();
            const edgeList& edgesA = cycPatch.edges();

            const labelList& meshPtsB = nbrPatch.meshPoints();
            const edgeList& edgesB = nbrPatch.edges();

            // Extract local values. Create map from edge to value.
            Map<T> half0Values(edgesA.size() / 20);
            Map<T> half1Values(half0Values.size());

            forAll(coupledEdges, edgei)
            {
                const edge& twoEdges = coupledEdges[edgei];

                {
                    const edge& e0 = edgesA[twoEdges[0]];
                    const edge meshEdge0(meshPtsA, e0);

                    const auto iter = edgeValues.cfind(meshEdge0);

                    if (iter.good())
                    {
                        half0Values.insert(edgei, iter.val());
                    }
                }
                {
                    const edge& e1 = edgesB[twoEdges[1]];
                    const edge meshEdge1(meshPtsB, e1);

                    const auto iter = edgeValues.cfind(meshEdge1);

                    if (iter.good())
                    {
                        half1Values.insert(edgei, iter.val());
                    }
                }
            }

            // Transform to this side
            top(cycPatch, half1Values);
            top(nbrPatch, half0Values);


            // Extract and combine information

            forAll(coupledEdges, edgei)
            {
                const edge& twoEdges = coupledEdges[edgei];

                const auto half1Fnd = half1Values.cfind(edgei);

                if (half1Fnd.good())
                {
                    const edge& e0 = edgesA[twoEdges[0]];
                    const edge meshEdge0(meshPtsA, e0);

                    combine
                    (
                        edgeValues,
                        cop,
                        meshEdge0,      // edge
                        *half1Fnd       // value
                    );
                }

                const auto half0Fnd = half0Values.cfind(edgei);

                if (half0Fnd.good())
                {
                    const edge& e1 = edgesB[twoEdges[1]];
                    const edge meshEdge1(meshPtsB, e1);

                    combine
                    (
                        edgeValues,
                        cop,
                        meshEdge1,      // edge
                        *half0Fnd       // value
                    );
                }
            }
        }
    }

    // Synchronize multiple shared points.
    // Problem is that we don't want to construct shared edges so basically
    // we do here like globalMeshData but then using sparse edge representation
    // (EdgeMap instead of mesh.edges())

    const globalMeshData& pd = mesh.globalData();
    const labelList& sharedPtAddr = pd.sharedPointAddr();
    const labelList& sharedPtLabels = pd.sharedPointLabels();

    // 1. Create map from meshPoint to globalShared index.
    Map<label> meshToShared(sharedPtLabels, sharedPtAddr);

    // Values on shared points. From two sharedPtAddr indices to a value.
    EdgeMap<T> sharedEdgeValues(meshToShared.size());

    // From shared edge to mesh edge. Used for merging later on.
    EdgeMap<edge> potentialSharedEdge(meshToShared.size());

    // 2. Find any edges using two global shared points. These will always be
    // on the outside of the mesh. (though might not be on coupled patch
    // if is single edge and not on coupled face)
    // Store value (if any) on sharedEdgeValues
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); ++facei)
    {
        const face& f = mesh.faces()[facei];

        forAll(f, fp)
        {
            const label v0 = f[fp];
            const label v1 = f[f.fcIndex(fp)];

            const auto v0Fnd = meshToShared.cfind(v0);

            if (v0Fnd.good())
            {
                const auto v1Fnd = meshToShared.cfind(v1);

                if (v1Fnd.good())
                {
                    const edge meshEdge(v0, v1);

                    // edge in shared point labels
                    const edge sharedEdge(*v0Fnd, *v1Fnd);

                    // Store mesh edge as a potential shared edge.
                    potentialSharedEdge.insert(sharedEdge, meshEdge);

                    const auto edgeFnd = edgeValues.cfind(meshEdge);

                    if (edgeFnd.good())
                    {
                        // edge exists in edgeValues. See if already in map
                        // (since on same processor, e.g. cyclic)
                        combine
                        (
                            sharedEdgeValues,
                            cop,
                            sharedEdge, // edge
                            *edgeFnd    // value
                        );
                    }
                }
            }
        }
    }


    // Now sharedEdgeValues will contain per potential sharedEdge the value.
    // (potential since an edge having two shared points is not necessary a
    //  shared edge).
    // Reduce this on the master.

    if (UPstream::parRun())
    {
        if (UPstream::master())
        {
            // Receive the edges using shared points from other procs
            for (const int proci : UPstream::subProcs())
            {
                EdgeMap<T> nbrValues;
                IPstream::recv(nbrValues, proci);

                // Merge neighbouring values with my values
                forAllConstIters(nbrValues, iter)
                {
                    combine
                    (
                        sharedEdgeValues,
                        cop,
                        iter.key(),     // edge
                        iter.val()      // value
                    );
                }
            }
        }
        else
        {
            // Send to master
            OPstream::send(sharedEdgeValues, UPstream::masterNo());
        }

        // Broadcast: send merged values to all
        Pstream::broadcast(sharedEdgeValues);
    }


    // Merge sharedEdgeValues (keyed on sharedPointAddr) into edgeValues
    // (keyed on mesh points).

    // Loop over all my shared edges.
    forAllConstIters(potentialSharedEdge, iter)
    {
        const edge& sharedEdge = iter.key();
        const edge& meshEdge = iter.val();

        // Do I have a value for the shared edge?
        const auto sharedFnd = sharedEdgeValues.cfind(sharedEdge);

        if (sharedFnd.good())
        {
            combine
            (
                edgeValues,
                cop,
                meshEdge,       // edge
                *sharedFnd      // value
            );
        }
    }

    // Reset tag
    UPstream::msgType(oldTag);
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    List<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " != number of points " << mesh.nPoints() << nl
            << abort(FatalError);
    }

    mesh.globalData().syncPointData(pointValues, cop, top);
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    const labelUList& meshPoints,
    List<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (pointValues.size() != meshPoints.size())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " != number of meshPoints " << meshPoints.size() << nl
            << abort(FatalError);
    }
    const globalMeshData& gd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gd.coupledPatch();
    const Map<label>& mpm = cpp.meshPointMap();

    List<T> cppFld(cpp.nPoints(), nullValue);

    forAll(meshPoints, i)
    {
        const auto iter = mpm.cfind(meshPoints[i]);

        if (iter.good())
        {
            cppFld[iter.val()] = pointValues[i];
        }
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalPointSlaves(),
        gd.globalPointTransformedSlaves(),
        gd.globalPointSlavesMap(),
        gd.globalTransforms(),
        cop,
        top
    );

    forAll(meshPoints, i)
    {
        const auto iter = mpm.cfind(meshPoints[i]);

        if (iter.good())
        {
            pointValues[i] = cppFld[iter.val()];
        }
    }
}


template<class T, class CombineOp, class TransformOp, class FlipOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    List<T>& edgeValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top,
    const FlipOp& fop
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " != number of edges " << mesh.nEdges() << nl
            << abort(FatalError);
    }

    const edgeList& edges = mesh.edges();
    const globalMeshData& gd = mesh.globalData();
    const labelList& meshEdges = gd.coupledPatchMeshEdges();
    const indirectPrimitivePatch& cpp = gd.coupledPatch();
    const edgeList& cppEdges = cpp.edges();
    const labelList& mp = cpp.meshPoints();
    const globalIndexAndTransform& git = gd.globalTransforms();
    const mapDistribute& edgeMap = gd.globalEdgeSlavesMap();
    const bitSet& orientation = gd.globalEdgeOrientation();

    List<T> cppFld(meshEdges.size());
    forAll(meshEdges, i)
    {
        const edge& cppE = cppEdges[i];
        const label meshEdgei = meshEdges[i];
        const edge& meshE = edges[meshEdgei];

        // 1. is cpp edge oriented as mesh edge
        // 2. is cpp edge oriented same as master edge

        const int dir = edge::compare(meshE, edge(mp, cppE));
        if (dir == 0)
        {
            FatalErrorInFunction<< "Problem:"
                << " mesh edge:" << meshE.line(mesh.points())
                << " coupled edge:" << cppE.line(cpp.localPoints())
                << exit(FatalError);
        }

        const bool sameOrientation = ((dir == 1) == orientation[i]);

        if (sameOrientation)
        {
            cppFld[i] = edgeValues[meshEdgei];
        }
        else
        {
            cppFld[i] = fop(edgeValues[meshEdgei]);
        }
    }


    globalMeshData::syncData
    (
        cppFld,
        gd.globalEdgeSlaves(),
        gd.globalEdgeTransformedSlaves(),
        edgeMap,
        git,
        cop,
        top
    );

    // Extract back onto mesh
    forAll(meshEdges, i)
    {
        const edge& cppE = cppEdges[i];
        const label meshEdgei = meshEdges[i];
        const edge& meshE = edges[meshEdgei];

        // 1. is cpp edge oriented as mesh edge
        // 2. is cpp edge oriented same as master edge

        const int dir = edge::compare(meshE, edge(mp, cppE));
        const bool sameOrientation = ((dir == 1) == orientation[i]);

        if (sameOrientation)
        {
            edgeValues[meshEdges[i]] = cppFld[i];
        }
        else
        {
            edgeValues[meshEdges[i]] =  fop(cppFld[i]);
        }
    }
}


template<class T, class CombineOp, class TransformOp, class FlipOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    const labelUList& meshEdges,
    List<T>& edgeValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top,
    const FlipOp& fop
)
{
    if (edgeValues.size() != meshEdges.size())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " != number of meshEdges " << meshEdges.size() << nl
            << abort(FatalError);
    }
    const edgeList& edges = mesh.edges();
    const globalMeshData& gd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gd.coupledPatch();
    const edgeList& cppEdges = cpp.edges();
    const labelList& mp = cpp.meshPoints();
    const Map<label>& mpm = gd.coupledPatchMeshEdgeMap();
    const bitSet& orientation = gd.globalEdgeOrientation();

    List<T> cppFld(cpp.nEdges(), nullValue);

    forAll(meshEdges, i)
    {
        const label meshEdgei = meshEdges[i];
        const auto iter = mpm.cfind(meshEdgei);
        if (iter.good())
        {
            const label cppEdgei = iter.val();
            const edge& cppE = cppEdges[cppEdgei];
            const edge& meshE = edges[meshEdgei];

            // 1. is cpp edge oriented as mesh edge
            // 2. is cpp edge oriented same as master edge

            const int dir = edge::compare(meshE, edge(mp, cppE));
            if (dir == 0)
            {
                FatalErrorInFunction<< "Problem:"
                    << " mesh edge:" << meshE.line(mesh.points())
                    << " coupled edge:" << cppE.line(cpp.localPoints())
                    << exit(FatalError);
            }

            const bool sameOrientation = ((dir == 1) == orientation[i]);

            if (sameOrientation)
            {
                cppFld[cppEdgei] = edgeValues[i];
            }
            else
            {
                cppFld[cppEdgei] = fop(edgeValues[i]);
            }
        }
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalEdgeSlaves(),
        gd.globalEdgeTransformedSlaves(),
        gd.globalEdgeSlavesMap(),
        gd.globalTransforms(),
        cop,
        top
    );

    forAll(meshEdges, i)
    {
        label meshEdgei = meshEdges[i];
        const auto iter = mpm.cfind(meshEdgei);
        if (iter.good())
        {
            label cppEdgei = iter.val();
            const edge& cppE = cppEdges[cppEdgei];
            const edge& meshE = edges[meshEdgei];

            const int dir = edge::compare(meshE, edge(mp, cppE));
            const bool sameOrientation = ((dir == 1) == orientation[i]);

            if (sameOrientation)
            {
                edgeValues[i] = cppFld[cppEdgei];
            }
            else
            {
                edgeValues[i] = fop(cppFld[cppEdgei]);
            }
        }
    }
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncBoundaryFaceList
(
    const polyMesh& mesh,
    UList<T>& faceValues,
    const CombineOp& cop,
    const TransformOp& top,
    const bool parRun
)
{
    // Offset (global to local) for start of boundaries
    const label boundaryOffset = mesh.nInternalFaces();

    if (faceValues.size() != mesh.nBoundaryFaces())
    {
        FatalErrorInFunction
            << "Number of values " << faceValues.size()
            << " != number of boundary faces " << mesh.nBoundaryFaces() << nl
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    if (parRun && UPstream::parRun())
    {
        // Avoid mesh.globalData() - possible race condition

        if
        (
            is_contiguous<T>::value
         && UPstream::defaultCommsType == UPstream::commsTypes::nonBlocking
        )
        {
            const label startRequest = UPstream::nRequests();

            // Receive buffer
            List<T> receivedValues(mesh.nBoundaryFaces());

            // Set up reads
            for (const polyPatch& pp : patches)
            {
                const auto* ppp = isA<processorPolyPatch>(pp);

                if (ppp && pp.size())
                {
                    const auto& procPatch = *ppp;

                    SubList<T> fld
                    (
                        receivedValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );

                    UIPstream::read
                    (
                        UPstream::commsTypes::nonBlocking,
                        procPatch.neighbProcNo(),
                        fld.data_bytes(),
                        fld.size_bytes()
                    );
                }
            }

            // Set up writes
            for (const polyPatch& pp : patches)
            {
                const auto* ppp = isA<processorPolyPatch>(pp);

                if (ppp && pp.size())
                {
                    const auto& procPatch = *ppp;

                    const SubList<T> fld
                    (
                        faceValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );

                    UOPstream::write
                    (
                        UPstream::commsTypes::nonBlocking,
                        procPatch.neighbProcNo(),
                        fld.cdata_bytes(),
                        fld.size_bytes()
                    );
                }
            }

            // Wait for all comms to finish
            UPstream::waitRequests(startRequest);

            // Combine with existing data
            for (const polyPatch& pp : patches)
            {
                const auto* ppp = isA<processorPolyPatch>(pp);

                if (ppp && pp.size())
                {
                    const auto& procPatch = *ppp;

                    SubList<T> recvFld
                    (
                        receivedValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );
                    const List<T>& fakeList = recvFld;
                    top(procPatch, const_cast<List<T>&>(fakeList));

                    SubList<T> patchValues
                    (
                        faceValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );

                    forAll(patchValues, i)
                    {
                        cop(patchValues[i], recvFld[i]);
                    }
                }
            }
        }
        else
        {
            DynamicList<label> neighbProcs;
            PstreamBuffers pBufs;

            // Send
            for (const polyPatch& pp : patches)
            {
                const auto* ppp = isA<processorPolyPatch>(pp);

                if (ppp && pp.size())
                {
                    const auto& procPatch = *ppp;
                    const label nbrProci = procPatch.neighbProcNo();

                    // Neighbour connectivity
                    neighbProcs.push_uniq(nbrProci);

                    const SubList<T> fld
                    (
                        faceValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );

                    UOPstream toNbr(nbrProci, pBufs);
                    toNbr << fld;
                }
            }

            // Limit exchange to involved procs
            pBufs.finishedNeighbourSends(neighbProcs);


            // Receive and combine.
            for (const polyPatch& pp : patches)
            {
                const auto* ppp = isA<processorPolyPatch>(pp);

                if (ppp && pp.size())
                {
                    const auto& procPatch = *ppp;
                    const label nbrProci = procPatch.neighbProcNo();

                    List<T> recvFld;
                    {
                        UIPstream fromNbr(nbrProci, pBufs);
                        fromNbr >> recvFld;
                    }

                    top(procPatch, recvFld);

                    SubList<T> patchValues
                    (
                        faceValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );

                    forAll(patchValues, i)
                    {
                        cop(patchValues[i], recvFld[i]);
                    }
                }
            }
        }
    }

    // Do the cyclics.
    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            // Owner does all.

            const cyclicPolyPatch& cycPatch = *cpp;
            const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
            const label patchSize = cycPatch.size();

            SubList<T> ownPatchValues
            (
                faceValues,
                patchSize,
                cycPatch.start()-boundaryOffset
            );

            SubList<T> nbrPatchValues
            (
                faceValues,
                patchSize,
                nbrPatch.start()-boundaryOffset
            );

            // Transform (copy of) data on both sides
            List<T> ownVals(ownPatchValues);
            top(nbrPatch, ownVals);

            List<T> nbrVals(nbrPatchValues);
            top(cycPatch, nbrVals);

            forAll(ownPatchValues, i)
            {
                cop(ownPatchValues[i], nbrVals[i]);
            }

            forAll(nbrPatchValues, i)
            {
                cop(nbrPatchValues[i], ownVals[i]);
            }
        }
    }

    // Reset tag
    UPstream::msgType(oldTag);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned Width, class CombineOp>
void Foam::syncTools::syncFaceList
(
    const polyMesh& mesh,
    const bool isBoundaryOnly,
    PackedList<Width>& faceValues,
    const CombineOp& cop,
    const bool parRun
)
{
    // Offset (global to local) for start of boundaries
    const label boundaryOffset = (isBoundaryOnly ? mesh.nInternalFaces() : 0);

    // Check size
    if (faceValues.size() != (mesh.nFaces() - boundaryOffset))
    {
        FatalErrorInFunction
            << "Number of values " << faceValues.size()
            << " != number of "
            << (isBoundaryOnly ? "boundary" : "mesh") << " faces "
            << (mesh.nFaces() - boundaryOffset) << nl
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    if (parRun && UPstream::parRun())
    {
        const label startRequest = UPstream::nRequests();

        // Receive buffers
        PtrList<PackedList<Width>> recvBufs(patches.size());

        // Set up reads
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.size())
            {
                const auto& procPatch = *ppp;
                const label patchi = pp.index();

                auto& recvbuf = recvBufs.emplace_set(patchi, pp.size());

                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    procPatch.neighbProcNo(),
                    recvbuf.data_bytes(),
                    recvbuf.size_bytes()
                );
            }
        }

        // Send buffers
        PtrList<PackedList<Width>> sendBufs(patches.size());

        // Set up writes
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.size())
            {
                const auto& procPatch = *ppp;
                const label patchi = pp.index();

                const labelRange range(pp.start()-boundaryOffset, pp.size());

                auto& sendbuf = sendBufs.emplace_set(patchi, faceValues, range);

                UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    procPatch.neighbProcNo(),
                    sendbuf.cdata_bytes(),
                    sendbuf.size_bytes()
                );
            }
        }

        // Wait for all comms to finish
        UPstream::waitRequests(startRequest);

        // Combine with existing data
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.size())
            {
                const label patchi = pp.index();
                const label patchSize = pp.size();

                const auto& recvbuf = recvBufs[patchi];

                // Combine (bitwise)
                label bFacei = pp.start()-boundaryOffset;
                for (label i = 0; i < patchSize; ++i)
                {
                    unsigned int recvVal = recvbuf[i];
                    unsigned int faceVal = faceValues[bFacei];

                    cop(faceVal, recvVal);
                    faceValues.set(bFacei, faceVal);

                    ++bFacei;
                }
            }
        }
    }


    // Do the cyclics.
    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            // Owner does all.

            const cyclicPolyPatch& cycPatch = *cpp;
            const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
            const label patchSize = cycPatch.size();

            label face0 = cycPatch.start()-boundaryOffset;
            label face1 = nbrPatch.start()-boundaryOffset;
            for (label i = 0; i < patchSize; ++i)
            {
                unsigned int val0 = faceValues[face0];
                unsigned int val1 = faceValues[face1];

                unsigned int t = val0;
                cop(t, val1);
                faceValues[face0] = t;

                cop(val1, val0);
                faceValues[face1] = val1;

                ++face0;
                ++face1;
            }
        }
    }

    // Reset tag
    UPstream::msgType(oldTag);
}


template<class T>
void Foam::syncTools::swapBoundaryCellList
(
    const polyMesh& mesh,
    const UList<T>& cellData,
    List<T>& neighbourCellData,
    const bool parRun
)
{
    if (cellData.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Number of cell values " << cellData.size()
            << " != number of cells " << mesh.nCells() << nl
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    neighbourCellData.resize(mesh.nBoundaryFaces());

    for (const polyPatch& pp : patches)
    {
        const auto& faceCells = pp.faceCells();

        // ie, boundarySlice() = patchInternalList()
        SubList<T>
        (
            neighbourCellData,
            faceCells.size(),
            pp.offset()
        ) = UIndirectList<T>(cellData, faceCells);
    }

    syncTools::swapBoundaryFaceList(mesh, neighbourCellData, parRun);
}


template<unsigned Width, class CombineOp>
void Foam::syncTools::syncFaceList
(
    const polyMesh& mesh,
    PackedList<Width>& faceValues,
    const CombineOp& cop,
    const bool parRun
)
{
    syncFaceList(mesh, false, faceValues, cop, parRun);
}


template<unsigned Width, class CombineOp>
void Foam::syncTools::syncBoundaryFaceList
(
    const polyMesh& mesh,
    PackedList<Width>& faceValues,
    const CombineOp& cop,
    const bool parRun
)
{
    syncFaceList(mesh, true, faceValues, cop, parRun);
}


template<unsigned Width>
void Foam::syncTools::swapFaceList
(
    const polyMesh& mesh,
    PackedList<Width>& faceValues,
    const bool parRun
)
{
    syncFaceList
    (
        mesh,
        faceValues,
        eqOp<unsigned int>(),
        parRun
    );
}


template<unsigned Width>
void Foam::syncTools::swapBoundaryFaceList
(
    const polyMesh& mesh,
    PackedList<Width>& faceValues,
    const bool parRun
)
{
    syncBoundaryFaceList
    (
        mesh,
        faceValues,
        eqOp<unsigned int>(),
        parRun
    );
}


template<unsigned Width, class CombineOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    PackedList<Width>& pointValues,
    const CombineOp& cop,
    const unsigned int nullValue
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " != number of points " << mesh.nPoints() << nl
            << abort(FatalError);
    }

    const globalMeshData& gd = mesh.globalData();
    const labelList& meshPoints = gd.coupledPatch().meshPoints();

    List<unsigned int> cppFld(gd.globalPointSlavesMap().constructSize());
    forAll(meshPoints, i)
    {
        cppFld[i] = pointValues[meshPoints[i]];
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalPointSlaves(),
        gd.globalPointTransformedSlaves(),
        gd.globalPointSlavesMap(),
        cop
    );

    // Extract back to mesh
    forAll(meshPoints, i)
    {
        pointValues[meshPoints[i]] = cppFld[i];
    }
}


template<unsigned Width, class CombineOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    PackedList<Width>& edgeValues,
    const CombineOp& cop,
    const unsigned int nullValue
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " != number of edges " << mesh.nEdges() << nl
            << abort(FatalError);
    }

    const globalMeshData& gd = mesh.globalData();
    const labelList& meshEdges = gd.coupledPatchMeshEdges();

    List<unsigned int> cppFld(gd.globalEdgeSlavesMap().constructSize());
    forAll(meshEdges, i)
    {
        cppFld[i] = edgeValues[meshEdges[i]];
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalEdgeSlaves(),
        gd.globalEdgeTransformedSlaves(),
        gd.globalEdgeSlavesMap(),
        cop
    );

    // Extract back to mesh
    forAll(meshEdges, i)
    {
        edgeValues[meshEdges[i]] = cppFld[i];
    }
}


// ************************************************************************* //
