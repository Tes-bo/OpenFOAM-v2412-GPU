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

#include "faMeshDistributor.H"
#include "globalIndex.H"
#include "BitOps.H"
#include "ListOps.H"
#include "mapDistributePolyMesh.H"
#include "processorFaPatch.H"
#include "labelPairHashes.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::mapDistributePolyMesh
Foam::faMeshDistributor::distribute
(
    const faMesh& oldMesh,
    const mapDistributePolyMesh& distMap,
    const polyMesh& tgtPolyMesh,
    autoPtr<faMesh>& newMeshPtr
)
{
    newMeshPtr.reset(nullptr);

    const uindirectPrimitivePatch& oldAreaPatch = oldMesh.patch();

    // Original (patch) sizes before any changes
    const label nOldPoints = oldAreaPatch.nPoints();
    const label nOldEdges = oldAreaPatch.nEdges();
    const label nOldFaces = oldAreaPatch.nFaces();


    // ------------------------
    // Step 1: The face mapping
    // ------------------------
    // Relatively straightforward.
    // - a subset of face selections on the polyMesh faceMap

    mapDistributeBase faFaceMap(distMap.faceMap());

    {
        // Retain faces needed for the faMesh - preserve original order

        // Note can use compactLocalData without problem since
        // finiteArea is only defined on real boundary faces so there
        // is no danger of sending an internal or processor face.

        labelList oldToNewSub;
        labelList oldToNewConstruct;

        faFaceMap.compactLocalData
        (
            oldMesh.faceLabels(),
            oldToNewSub,
            oldToNewConstruct,
            distMap.nOldFaces(),
            UPstream::msgType()
        );


        // The receiving side:
        // Mapped faces (>= 0) are the polyMesh face labels that define
        // the faMesh. Since the compact order is implicitly sorted in
        // ascending order, it tends to form contiguous ranges on the
        // polyPatches, which serves our purpose nicely.

        labelList newFaceLabels
        (
            ListOps::findIndices(oldToNewConstruct, labelRange::ge0())
        );

        // Set up to read-if-present
        IOobject io(tgtPolyMesh);
        io.readOpt(IOobject::READ_IF_PRESENT);

        newMeshPtr.reset
        (
            new faMesh(tgtPolyMesh, std::move(newFaceLabels), io)
        );
    }

    // The face map is now complete.

    // The new faMesh and the corresponding primitive patch
    auto& newMesh = newMeshPtr();
    const uindirectPrimitivePatch& newAreaPatch = newMesh.patch();


    // ------------------------
    // Step 2: The edge mapping
    // ------------------------
    // Use globally unique addressing to identify both sides of the edges.

    // A local bookkeeping struct of (face0 face1 edge0 edge1)
    // selected for a stable and unique sort

    struct faceEdgeTuple : public FixedList<label, 4>
    {
        // Inherit constructors
        using FixedList<label, 4>::FixedList;

        // Default construct as 'invalid'
        faceEdgeTuple() : FixedList<label, 4>(-1) {}

        // Is empty if face indices are negative
        bool empty() const
        {
            return (face0() < 0 && face1() < 0);
        }

        // Global face numbers are the first sort index
        label face0() const { return (*this)[0]; }
        label face1() const { return (*this)[1]; }

        // Additional sorting based on edges
        label edge0() const { return (*this)[2]; }
        label edge1() const { return (*this)[3]; }

        label getFace(int i) const { return (*this)[(i ? 1 : 0)]; }
        label getEdge(int i) const { return (*this)[(i ? 3 : 2)]; }

        labelPair getFaces() const { return labelPair(face0(), face1()); }
        labelPair getPair0() const { return labelPair(face0(), edge0()); }
        // labelPair getPair1() const { return labelPair(face1(), edge1()); }

        // Canonically sorted order
        void canonical()
        {
            if (face1() >= 0 && face0() >= face1())
            {
                std::swap((*this)[0], (*this)[1]);
                std::swap((*this)[2], (*this)[3]);
            }
        }

        // Slot a face-edge pair into a free location in the tuple
        void add(const labelPair& faceEdge)
        {
            if ((*this)[0] < 0)  // face0
            {
                (*this)[0] = faceEdge.first();  // face
                (*this)[2] = faceEdge.second(); // edge
            }
            else if ((*this)[1] < 0)  // face1
            {
                (*this)[1] = faceEdge.first();  // face
                (*this)[3] = faceEdge.second(); // edge
            }
        }

        // Combine operation
        void combine(const faceEdgeTuple& y)
        {
            auto& x = *this;

            if (y.empty() || x == y)
            {
                // Nothing to do
            }
            else if (x.empty())
            {
                x = y;
            }
            else  // Both non-empty, but non-identical (should not occur)
            {
                FatalErrorInFunction
                    << "Unexpected edge matching: "
                    << x << " vs. " << y << endl
                    << exit(FatalError);
            }
        }

        // Combine operation
        void operator()(faceEdgeTuple& x, const faceEdgeTuple& y) const
        {
            x.combine(y);
        }
    };


    // ------------------------
    // (Edge mapping)
    // ------------------------
    // 1.
    // Create a list of destination edges for each face,
    // appended by a unique face identifier.
    // Using global face numbering from the *target* mesh.

    const globalIndex uniqFaceIndex(newAreaPatch.nFaces());

    labelListList dstFaceEdgeIds(newAreaPatch.nFaces());

    forAll(newAreaPatch.faceEdges(), facei)
    {
        const labelList& fEdges = newAreaPatch.faceEdges()[facei];
        labelList& dstEdges = dstFaceEdgeIds[facei];

        dstEdges.resize(fEdges.size() + 1);

        // Leading entries are the destination (patch) edge indices
        labelSubList(dstEdges, fEdges.size()) = fEdges;

        // Last entry is the face id
        dstEdges.last() = uniqFaceIndex.toGlobal(facei);
    }

    // Send back to the original mesh locations
    faFaceMap.reverseDistribute(nOldFaces, dstFaceEdgeIds);


    // 2.
    // Walk all original faces and their edges to generate a edge lookup
    // table with the destination face/edge information.
    // Eg ((globFacei, localEdgei), (globFacei, localEdgei))

    // NB: currently no provision for face flips, which would potentially
    // reverse the local edge order.

    List<faceEdgeTuple> halfEdgeLookup(nOldEdges, faceEdgeTuple());

    forAll(oldAreaPatch.faceEdges(), facei)
    {
        const labelList& fEdges = oldAreaPatch.faceEdges()[facei];

        // The corresponding destination edges (+ uniqFacei).
        const labelList& dstFaceEdges = dstFaceEdgeIds[facei];

        // Last entry has the unique destination face id
        const label uniqFacei = dstFaceEdges.last();

        forAll(fEdges, faceEdgei)
        {
            const label srcEdgei = fEdges[faceEdgei];
            const label dstEdgei = dstFaceEdges[faceEdgei];

            halfEdgeLookup[srcEdgei].add(labelPair(uniqFacei, dstEdgei));
        }
    }

    // 3.
    // Add patch indices - encoded as '-(patchId + 2)' for
    // non-processor boundaries (will be needed later).

    // Also record which procs get patches from here [proc] -> [patches]
    // Use for building patchMap
    Map<labelHashSet> patchMapInfo;

    label nNonProcessor(0);
    {
        const faBoundaryMesh& oldBndMesh = oldMesh.boundary();

        forAll(oldBndMesh, patchi)
        {
            const faPatch& fap = oldBndMesh[patchi];
            if (isA<processorFaPatch>(fap))
            {
                break;
            }
            ++nNonProcessor;

            // Encoded as -(patchId + 2)
            const labelPair encodedPatchId(-(patchi+2), -(patchi+2));

            for (const label srcEdgei : fap.edgeLabels())
            {
                faceEdgeTuple& facePairings = halfEdgeLookup[srcEdgei];

                facePairings.add(encodedPatchId);

                label dstFacei = facePairings.face0();

                label proci = uniqFaceIndex.whichProcID(dstFacei);

                patchMapInfo(proci).insert(patchi);
            }
        }

        Pstream::broadcast(nNonProcessor);
    }

    // Processor-processor connections
    if (Pstream::parRun())
    {
        const label startOfRequests = UPstream::nRequests();

        const faBoundaryMesh& oldBndMesh = oldMesh.boundary();

        // Setup sends
        for (const faPatch& fap : oldBndMesh)
        {
            const auto* fapp = isA<processorFaPatch>(fap);

            if (fapp)
            {
                // Send (dstFacei dstEdgei) per patch edge.
                // Since we are walking a boundary edge, the first location
                // from the previously established lookup provides
                // our information.

                List<labelPair> edgeTuples(fap.edgeLabels().size());

                label edgeTuplei = 0;

                for (const label edgei : fap.edgeLabels())
                {
                    edgeTuples[edgeTuplei] = halfEdgeLookup[edgei].getPair0();
                    ++edgeTuplei;
                }

                fapp->send<labelPair>
                (
                    Pstream::commsTypes::nonBlocking,
                    edgeTuples
                );
            }
        }

        // Wait for outstanding requests
        // (commsType == UPstream::commsTypes::nonBlocking)
        UPstream::waitRequests(startOfRequests);

        // Receive values
        for (const faPatch& fap : oldBndMesh)
        {
            const auto* fapp = isA<processorFaPatch>(fap);

            if (fapp)
            {
                // Receive (dstFacei dstEdgei) per patch edge.
                // - slot into the second location of the lookup

                List<labelPair> edgeTuples(fap.edgeLabels().size());

                fapp->receive<labelPair>
                (
                    Pstream::commsTypes::nonBlocking,
                    edgeTuples
                );

                label edgeTuplei = 0;

                for (const label srcEdgei : fap.edgeLabels())
                {
                    halfEdgeLookup[srcEdgei].add(edgeTuples[edgeTuplei]);
                    ++edgeTuplei;
                }
            }
        }
    }

    // Globally consistent order
    for (auto& tuple : halfEdgeLookup)
    {
        tuple.canonical();
    }


    // Now ready to actually construct the edgeMap

    mapDistributeBase faEdgeMap
    (
        newAreaPatch.nEdges(),  // constructSize
        labelListList(),        // subMap
        labelListList(),        // constructMap
        false,                  // subHasFlip
        false,                  // constructHasFlip
        faFaceMap.comm()
    );

    {
        // Pass 1.
        // Count the number of edges to be sent to each proc

        labelList nProcEdges(Pstream::nProcs(faFaceMap.comm()), Zero);

        forAll(halfEdgeLookup, srcEdgei)
        {
            const faceEdgeTuple& facePairings = halfEdgeLookup[srcEdgei];

            labelPair dstProcs(-1, -1);

            for (int sidei = 0; sidei < 2; ++sidei)
            {
                const label dstFacei = facePairings.getFace(sidei);
                //const label dstEdgei = facePairings.getEdge(sidei);

                if (dstFacei < 0)
                {
                    continue;
                }

                label proci = uniqFaceIndex.whichProcID(dstFacei);
                dstProcs[sidei] = proci;

                if (proci != -1 && dstProcs[0] != dstProcs[1])
                {
                    ++nProcEdges[proci];
                }
            }
        }

        auto& edgeSubMap = faEdgeMap.subMap();
        auto& edgeCnstrMap = faEdgeMap.constructMap();

        edgeSubMap.resize(nProcEdges.size());
        edgeCnstrMap.resize(nProcEdges.size());

        labelListList remoteEdges(nProcEdges.size());

        forAll(nProcEdges, proci)
        {
            edgeSubMap[proci].resize(nProcEdges[proci], -1);
            remoteEdges[proci].resize(nProcEdges[proci], -1);
        }


        // Pass 2.
        // Fill in the maps

        nProcEdges = Zero;  // Reset counter

        forAll(halfEdgeLookup, srcEdgei)
        {
            const faceEdgeTuple& facePairings = halfEdgeLookup[srcEdgei];

            labelPair dstProcs(-1, -1);

            for (int sidei = 0; sidei < 2; ++sidei)
            {
                const label dstFacei = facePairings.getFace(sidei);
                const label dstEdgei = facePairings.getEdge(sidei);

                if (dstFacei < 0)
                {
                    continue;
                }

                label proci = uniqFaceIndex.whichProcID(dstFacei);
                dstProcs[sidei] = proci;

                if (proci != -1 && dstProcs[0] != dstProcs[1])
                {
                    edgeSubMap[proci][nProcEdges[proci]] = srcEdgei;
                    remoteEdges[proci][nProcEdges[proci]] = dstEdgei;
                    ++nProcEdges[proci];
                }
            }
        }

        // The remoteEdges are what we know locally about what will be
        // received, but not what is actually received.
        // So need an all-to-all exchange

        Pstream::exchange<labelList, label>
        (
            remoteEdges,
            edgeCnstrMap,
            UPstream::msgType(),
            faEdgeMap.comm()
        );
    }

    // The edge map is now complete [in PrimitivePatch edge order]


    if (oldMesh.hasInternalEdgeLabels())
    {
        // If there are gaps in the edge numbering or the
        // internal edge labels are out of sequence would
        // have to use compactLocalData etc before sending
        // But just issue an error for now

        FatalErrorInFunction
            << "Originating faMesh has gaps in the edge addressing"
            << " this is currently unsupported"
            << abort(FatalError);
    }


    // ------------------------
    // Patch edge labels
    // ------------------------

    faPatchList newFaPatches;

    // Distribute the edge lookups.
    // Needs full version for the combine operation

    mapDistributeBase::distribute<faceEdgeTuple, faceEdgeTuple, identityOp>
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>::null(),
        faEdgeMap.constructSize(),

        faEdgeMap.subMap(),
        faEdgeMap.subHasFlip(),

        faEdgeMap.constructMap(),
        faEdgeMap.constructHasFlip(),

        halfEdgeLookup,
        faceEdgeTuple(),       // nullValue
        faceEdgeTuple(),       // CombineOp
        identityOp(),          // NegateOp
        UPstream::msgType(),
        faEdgeMap.comm()
    );

    {
        // Pass 1.
        // Count the number of edges for each patch type

        Map<label> nNonProcEdges(2*nNonProcessor);
        LabelPairMap<label> nProcEdges(2*nNonProcessor);

        forAll(halfEdgeLookup, edgei)
        {
            const auto& tuple = halfEdgeLookup[edgei];

            labelPair target(tuple.getFaces());

            if (target[1] < -1)
            {
                // Neighbour face was patchId encoded value
                label patchi = mag(target[1])-2;
                ++nNonProcEdges(patchi);
            }
            else if (target[0] >= 0 && target[1] >= 0)
            {
                // From global face to proc id
                target[0] = uniqFaceIndex.whichProcID(target[0]);
                target[1] = uniqFaceIndex.whichProcID(target[1]);

                // A connection between different procs but involving
                // myProc?
                if
                (
                    target[0] != target[1]
                 &&
                    (
                        target[0] == Pstream::myProcNo()
                     || target[1] == Pstream::myProcNo()
                    )
                )
                {
                    ++nProcEdges(target);
                }
            }
        }

        label nPatches = (nNonProcessor + nProcEdges.size());

        newFaPatches.resize(nPatches);

        labelList nEdgeLabels(nPatches, Zero);
        labelListList newEdgeLabels(nPatches);

        LabelPairMap<label> procPairToPatchId(2*nProcEdges.size());

        // Map processor-pair to index in patches
        {
            nPatches = nNonProcessor;

            for (const labelPair& twoProcs : nProcEdges.sortedToc())
            {
                procPairToPatchId.set(twoProcs, nPatches);
                ++nPatches;
            }
        }

        // Presize edgeLabels arrays
        {
            nPatches = 0;

            for (label patchi = 0; patchi < nNonProcessor; ++patchi)
            {
                label nLabels = nNonProcEdges.lookup(nPatches, 0);

                nEdgeLabels[nPatches] = nLabels;

                newEdgeLabels[nPatches].resize(nLabels);

                ++nPatches;
            }

            // Processor patches
            for (const labelPair& twoProcs : nProcEdges.sortedToc())
            {
                label nLabels = nProcEdges.lookup(twoProcs, 0);

                nEdgeLabels[nPatches] = nLabels;

                newEdgeLabels[nPatches].resize(nLabels);

                ++nPatches;
            }
        }

        nEdgeLabels = Zero;

        // Populate edgeLabels arrays - walk in canonically sorted
        // order to ensure that both sides of processor edges
        // correspond.

        const labelList order(Foam::sortedOrder(halfEdgeLookup));

        for (const label edgei : order)
        {
            const auto& tuple = halfEdgeLookup[edgei];

            labelPair target(tuple.getFaces());

            label patchi = -1;

            if (target[1] < -1)
            {
                // Neighbour face was patchId encoded value
                patchi = mag(target[1])-2;
            }
            else if (target[0] >= 0 && target[1] >= 0)
            {
                // From global face to proc id
                target[0] = uniqFaceIndex.whichProcID(target[0]);
                target[1] = uniqFaceIndex.whichProcID(target[1]);

                // A connection between different procs but involving
                // myProc?
                if
                (
                    target[0] != target[1]
                 &&
                    (
                        target[0] == Pstream::myProcNo()
                     || target[1] == Pstream::myProcNo()
                    )
                )
                {
                    patchi = procPairToPatchId.lookup(target, -1);
                }
            }

            if (patchi >= 0)
            {
                const label idx = nEdgeLabels[patchi]++;
                newEdgeLabels[patchi][idx] = edgei;
            }
        }


        // Clone all non-processor patches

        nPatches = 0;
        for (label patchi = 0; patchi < nNonProcessor; ++patchi)
        {
            newFaPatches.set
            (
                nPatches,
                oldMesh.boundary()[patchi].clone
                (
                    newMesh.boundary(),
                    newEdgeLabels[nPatches],  // edgeLabels
                    nPatches,
                    oldMesh.boundary()[patchi].ngbPolyPatchIndex()
                )
            );
            ++nPatches;
        }

        // Any processor patches
        for (const labelPair& twoProcs : nProcEdges.sortedToc())
        {
            label nbrProcNo =
            (
                twoProcs[1] == Pstream::myProcNo()
              ? twoProcs[0] : twoProcs[1]
            );

            newFaPatches.set
            (
                nPatches,
                new processorFaPatch
                (
                    newEdgeLabels[nPatches],  // edgeLabels
                    nPatches,
                    newMesh.boundary(),
                    -1,                       // nbrPolyPatchi
                    Pstream::myProcNo(),
                    nbrProcNo
                )
            );
            ++nPatches;
        }
    }


    newMesh.addFaPatches(newFaPatches);
    newMesh.init(true);


    // At this stage we now have a complete mapping overview in
    // terms of the PrimitivePatch edge ordering which now must be
    // changed into faMesh edge ordering.

    {
        labelList oldToNewSub;
        labelList oldToNewConstruct;

        // Assume we use all patch edges for the faMesh too.
        oldToNewSub.resize(oldAreaPatch.nEdges(), -1);
        oldToNewConstruct.resize(newAreaPatch.nEdges(), -1);

        // Remap old edges
        label nEdges = 0;

        // Internal edgeLabels
        for
        (
            label edgei = 0;
            edgei < oldAreaPatch.nInternalEdges();
            ++edgei
        )
        {
            oldToNewSub[edgei] = nEdges++;
        }

        // Boundary edgeLabels
        for (const faPatch& fap : oldMesh.boundary())
        {
            for (const label edgei : fap.edgeLabels())
            {
                oldToNewSub[edgei] = nEdges++;
            }
        }

        // New edges
        nEdges = 0;

        // Internal edgeLabels
        for
        (
            label edgei = 0;
            edgei < newAreaPatch.nInternalEdges();
            ++edgei
        )
        {
            oldToNewConstruct[edgei] = nEdges++;
        }

        // Boundary edgeLabels
        for (const faPatch& fap : newMesh.boundary())
        {
            for (const label edgei : fap.edgeLabels())
            {
                oldToNewConstruct[edgei] = nEdges++;
            }
        }

        mapDistributeBase::renumberMap
        (
            faEdgeMap.subMap(),
            oldToNewSub,
            faEdgeMap.subHasFlip()
        );

        mapDistributeBase::renumberMap
        (
            faEdgeMap.constructMap(),
            oldToNewConstruct,
            faEdgeMap.constructHasFlip()
        );
    }


    // ------------------------
    // Patch mapping
    // ------------------------

    mapDistributeBase faPatchMap
    (
        newMesh.boundary().size(),  // constructSize
        labelListList(),        // subMap
        labelListList(),        // constructMap
        false,                  // subHasFlip
        false,                  // constructHasFlip
        faFaceMap.comm()
    );

    // For patch maps, would normally transcribe from patchMapInfo
    // gathered earlier. However, existing practice (APR-2022) for
    // faMesh decomposition is to map all non-processor patches

    {
        // Map all non-processor patches
        const label nProcs = UPstream::nProcs(faPatchMap.comm());

        faPatchMap.subMap().resize(nProcs, identity(nNonProcessor));
        faPatchMap.constructMap().resize(nProcs, identity(nNonProcessor));
    }

    /// {
    ///     // Transcribe from patchMapInfo gathered earlier.
    ///     // - transform Map of labelHashSet to labelListList
    ///
    ///     labelListList sendToRemote(Pstream::nProcs(faPatchMap.comm()));
    ///
    ///     forAllConstIters(patchMapInfo, iter)
    ///     {
    ///         const label proci = iter.key();
    ///         sendToRemote[proci] = iter.val().sortedToc();
    ///     }
    ///
    ///     auto& patchSubMap = faPatchMap.subMap();
    ///     auto& patchCnstrMap = faPatchMap.constructMap();
    ///
    ///     patchSubMap = sendToRemote;
    ///     patchCnstrMap.resize(patchSubMap.size());
    ///
    ///     // Change sendToRemote into recv-from-remote by using
    ///     // all-to-all exchange
    ///
    ///     Pstream::exchange<labelList, label>
    ///     (
    ///         sendToRemote,
    ///         patchCnstrMap,
    ///         UPstream::msgType(),
    ///         faPatchMap.comm()
    ///     );
    /// }


    // ------------------------
    // Point mapping
    // ------------------------

    mapDistributeBase faPointMap(distMap.pointMap());

    {
        // Retain meshPoints needed for the faMesh - preserve original order
        // Need both sides (local/remote) for correct compaction maps
        // without dangling points.

        labelList oldToNewSub;
        labelList oldToNewConstruct;

        faPointMap.compactData
        (
            oldAreaPatch.meshPoints(),
            newAreaPatch.meshPoints(),
            oldToNewSub,
            oldToNewConstruct,
            distMap.nOldPoints(),
            UPstream::msgType()
        );
    }


    return mapDistributePolyMesh
    (
        // Mesh before changes
        nOldPoints,
        nOldEdges,          // area: nOldEdges (volume: nOldFaces)
        nOldFaces,          // area: nOldFaces (volume: nOldCells)

        labelList(oldMesh.boundary().patchStarts()),
        labelList(),        // oldPatchNMeshPoints [unused]

        mapDistribute(std::move(faPointMap)),
        mapDistribute(std::move(faEdgeMap)), // area: edgeMap (volume: faceMap)
        mapDistribute(std::move(faFaceMap)), // area: faceMap (volume: cellMap)
        mapDistribute(std::move(faPatchMap))
    );
}


Foam::mapDistributePolyMesh
Foam::faMeshDistributor::distribute
(
    const faMesh& oldMesh,
    const mapDistributePolyMesh& distMap,
    autoPtr<faMesh>& newMeshPtr
)
{
    return faMeshDistributor::distribute
    (
        oldMesh,
        distMap,
        oldMesh.mesh(),  // polyMesh
        newMeshPtr
    );
}


// ************************************************************************* //
