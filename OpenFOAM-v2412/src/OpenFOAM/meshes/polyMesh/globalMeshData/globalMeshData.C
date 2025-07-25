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

#include "globalMeshData.H"
#include "globalPoints.H"
#include "polyMesh.H"
#include "mapDistribute.H"
#include "labelIOList.H"
#include "mergePoints.H"
#include "processorPolyPatch.H"
#include "processorTopologyNew.H"
#include "globalIndexAndTransform.H"
#include "ListOps.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(globalMeshData, 0);

const scalar globalMeshData::matchTol_ = 1e-8;

template<>
class minEqOp<labelPair>
{
public:
    void operator()(labelPair& x, const labelPair& y) const
    {
        x[0] = min(x[0], y[0]);
        x[1] = min(x[1], y[1]);
    }
};
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::globalMeshData::initProcAddr()
{
    processorPatchIndices_.resize_nocopy(mesh_.boundaryMesh().size());
    processorPatchIndices_ = -1;

    processorPatchNeighbours_.resize_nocopy(mesh_.boundaryMesh().size());
    processorPatchNeighbours_ = -1;

    // Construct processor patch indexing. processorPatchNeighbours_ only
    // set if running in parallel!
    processorPatches_.resize_nocopy(mesh_.boundaryMesh().size());

    label nNeighbours = 0;

    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
        {
            processorPatches_[nNeighbours] = patchi;
            processorPatchIndices_[patchi] = nNeighbours++;
        }
    }
    processorPatches_.resize(nNeighbours);


    if (UPstream::parRun())
    {
        // Allocate unique tag for all comms
        const int oldTag = UPstream::incrMsgType();

        PstreamBuffers pBufs(mesh_.comm());

        // Send indices of my processor patches to my neighbours
        for (const label patchi : processorPatches_)
        {
            UOPstream toNeighbour
            (
                refCast<const processorPolyPatch>
                (
                    mesh_.boundaryMesh()[patchi]
                ).neighbProcNo(),
                pBufs
            );

            toNeighbour << processorPatchIndices_[patchi];
        }

        pBufs.finishedSends();

        for (const label patchi : processorPatches_)
        {
            UIPstream fromNeighbour
            (
                refCast<const processorPolyPatch>
                (
                    mesh_.boundaryMesh()[patchi]
                ).neighbProcNo(),
                pBufs
            );

            fromNeighbour >> processorPatchNeighbours_[patchi];
        }

        // Reset tag
        UPstream::msgType(oldTag);
    }
}


void Foam::globalMeshData::calcSharedPoints() const
{
    if
    (
        nGlobalPoints_ != -1
     || sharedPointLabelsPtr_
     || sharedPointAddrPtr_
    )
    {
        FatalErrorInFunction
            << "Shared point addressing already done" << abort(FatalError);
    }

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // Calculate all shared points (exclude points that are only
    // on two coupled patches). This does all the hard work.
    const globalPoints parallelPoints(mesh_, false, true);

    // Count the number of master points
    label nMaster = 0;
    forAll(parallelPoints.pointPoints(), i)
    {
        const labelList& pPoints = parallelPoints.pointPoints()[i];
        const labelList& transPPoints =
            parallelPoints.transformedPointPoints()[i];

        if (pPoints.size()+transPPoints.size() > 0)
        {
            nMaster++;
        }
    }

    // Allocate global numbers
    const globalIndex masterNumbering(nMaster);

    nGlobalPoints_ = masterNumbering.totalSize();


    // Push master number to slaves
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1. Fill master and slave slots
    nMaster = 0;
    labelList master(parallelPoints.map().constructSize(), -1);
    forAll(parallelPoints.pointPoints(), i)
    {
        const labelList& pPoints = parallelPoints.pointPoints()[i];
        const labelList& transPPoints =
            parallelPoints.transformedPointPoints()[i];

        if (pPoints.size()+transPPoints.size() > 0)
        {
            master[i] = masterNumbering.toGlobal(nMaster);

            labelUIndList(master, pPoints) = master[i];
            labelUIndList(master, transPPoints) = master[i];

            ++nMaster;
        }
    }


    // 2. Push slave slots back to local storage on originating processor
    // For all the four types of points:
    // - local master : already set
    // - local transformed slave point : the reverse transform at
    //   reverseDistribute will have copied it back to its originating local
    //   point
    // - remote untransformed slave point : sent back to originating processor
    // - remote transformed slave point : the reverse transform will
    //   copy it back into the remote slot which then gets sent back to
    //   originating processor

    parallelPoints.map().reverseDistribute
    (
        parallelPoints.map().constructSize(),
        master
    );


    // Collect all points that are a master or refer to a master.
    nMaster = 0;
    forAll(parallelPoints.pointPoints(), i)
    {
        if (master[i] != -1)
        {
            nMaster++;
        }
    }

    sharedPointLabelsPtr_.reset(new labelList(nMaster));
    labelList& sharedPointLabels = sharedPointLabelsPtr_();
    sharedPointAddrPtr_.reset(new labelList(nMaster));
    labelList& sharedPointAddr = sharedPointAddrPtr_();
    nMaster = 0;

    forAll(parallelPoints.pointPoints(), i)
    {
        if (master[i] != -1)
        {
            // I am master or slave
            sharedPointLabels[nMaster] = i;
            sharedPointAddr[nMaster] = master[i];
            nMaster++;
        }
    }

    // Reset tag
    UPstream::msgType(oldTag);

    if (debug)
    {
        Pout<< "globalMeshData : nGlobalPoints_:" << nGlobalPoints_ << nl
            << "globalMeshData : sharedPointLabels_:"
            << sharedPointLabelsPtr_().size() << nl
            << "globalMeshData : sharedPointAddr_:"
            << sharedPointAddrPtr_().size() << endl;
    }
}


void Foam::globalMeshData::countSharedEdges
(
    const EdgeMap<labelList>& procSharedEdges,
    EdgeMap<label>& globalShared,
    label& sharedEdgeI
)
{
    // Count occurrences of procSharedEdges in global shared edges table.
    forAllConstIters(procSharedEdges, iter)
    {
        const edge& e = iter.key();

        auto globalFnd = globalShared.find(e);

        if (globalFnd.good())
        {
            if (globalFnd() == -1)
            {
                // Second time occurrence of this edge.
                // Assign proper edge label.
                globalFnd() = sharedEdgeI++;
            }
        }
        else
        {
            // First time occurrence of this edge. Check how many we are adding.
            if (iter().size() == 1)
            {
                // Only one edge. Mark with special value.
                globalShared.insert(e, -1);
            }
            else
            {
                // Edge used more than once (even by local shared edges alone)
                // so allocate proper shared edge label.
                globalShared.insert(e, sharedEdgeI++);
            }
        }
    }
}


void Foam::globalMeshData::calcSharedEdges() const
{
    // Shared edges are shared between multiple processors. By their nature both
    // of their endpoints are shared points. (but not all edges using two shared
    // points are shared edges! There might e.g. be an edge between two
    // unrelated clusters of shared points)

    if
    (
        nGlobalEdges_ != -1
     || sharedEdgeLabelsPtr_
     || sharedEdgeAddrPtr_
    )
    {
        FatalErrorInFunction
            << "Shared edge addressing already done" << abort(FatalError);
    }


    const labelList& sharedPtAddr = sharedPointAddr();
    const labelList& sharedPtLabels = sharedPointLabels();

    // Since don't want to construct pointEdges for whole mesh create
    // Map for all shared points.
    Map<label> meshToShared(invertToMap(sharedPtLabels));


    // Find edges using shared points. Store correspondence to local edge
    // numbering. Note that multiple local edges can have the same shared
    // points! (for cyclics or separated processor patches)
    EdgeMap<labelList> localShared(2*sharedPtAddr.size());

    const edgeList& edges = mesh_.edges();

    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];

        const auto e0Fnd = meshToShared.cfind(e[0]);

        if (e0Fnd.good())
        {
            const auto e1Fnd = meshToShared.cfind(e[1]);

            if (e1Fnd.good())
            {
                // Found edge which uses shared points. Probably shared.

                // Construct the edge in shared points (or rather global indices
                // of the shared points)
                edge sharedEdge
                (
                    sharedPtAddr[e0Fnd.val()],
                    sharedPtAddr[e1Fnd.val()]
                );

                // Add this edge to list of edge labels
                localShared(sharedEdge).push_back(edgei);
            }
        }
    }


    // Now we have a table on every processors which gives its edges which use
    // shared points. Send this all to the master and have it allocate
    // global edge numbers for it. But only allocate a global edge number for
    // edge if it is used more than once!
    // Note that we are now sending the whole localShared to the master whereas
    // we only need the local count (i.e. the number of times a global edge is
    // used). But then this only gets done once so not too bothered about the
    // extra global communication.

    EdgeMap<label> globalShared(2*nGlobalPoints());

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    if (UPstream::master())
    {
        label sharedEdgeI = 0;

        // Merge my shared edges into the global list
        if (debug)
        {
            Pout<< "globalMeshData::calcSharedEdges : Merging in from proc0 : "
                << localShared.size() << endl;
        }
        countSharedEdges(localShared, globalShared, sharedEdgeI);

        // Receive data and insert
        if (UPstream::parRun())
        {
            for (const int proci : UPstream::subProcs())
            {
                // Receive the edges using shared points from the slave.
                EdgeMap<labelList> procSharedEdges;
                IPstream::recv(procSharedEdges, proci);

                if (debug)
                {
                    Pout<< "globalMeshData::calcSharedEdges : "
                        << "Merging in from proc"
                        << proci << " : " << procSharedEdges.size()
                        << endl;
                }
                countSharedEdges(procSharedEdges, globalShared, sharedEdgeI);
            }
        }

        // Now our globalShared should have some edges with -1 as edge label
        // These were only used once so are not proper shared edges.
        // Remove them.
        {
            EdgeMap<label> oldSharedEdges(std::move(globalShared));
            globalShared.clear();

            forAllConstIters(oldSharedEdges, iter)
            {
                if (iter.val() != -1)
                {
                    globalShared.insert(iter.key(), iter.val());
                }
            }

            if (debug)
            {
                Pout<< "globalMeshData::calcSharedEdges : Filtered "
                    << oldSharedEdges.size()
                    << " down to " << globalShared.size() << endl;
            }
        }
    }
    else
    {
        if (UPstream::parRun())
        {
            // send local edges to master
            OPstream::send(localShared, UPstream::masterNo());
        }
    }

    // Broadcast: merged edges to all
    Pstream::broadcast(globalShared);  // == worldComm;


    // Now use the global shared edges list (globalShared) to classify my local
    // ones (localShared)

    nGlobalEdges_ = globalShared.size();

    DynamicList<label> dynSharedEdgeLabels(globalShared.size());
    DynamicList<label> dynSharedEdgeAddr(globalShared.size());

    forAllConstIters(localShared, iter)
    {
        const edge& e = iter.key();

        const auto edgeFnd = globalShared.cfind(e);

        if (edgeFnd.good())
        {
            // My local edge is indeed a shared one. Go through all local edge
            // labels with this point combination.
            const labelList& edgeLabels = iter.val();

            for (const label edgei : edgeLabels)
            {
                // Store label of local mesh edge
                dynSharedEdgeLabels.append(edgei);

                // Store label of shared edge
                dynSharedEdgeAddr.append(edgeFnd());
            }
        }
    }


    sharedEdgeLabelsPtr_.reset
    (
        new labelList(std::move(dynSharedEdgeLabels))
    );

    sharedEdgeAddrPtr_.reset
    (
        new labelList(std::move(dynSharedEdgeAddr))
    );

    // Reset tag
    UPstream::msgType(oldTag);

    if (debug)
    {
        Pout<< "globalMeshData : nGlobalEdges_:" << nGlobalEdges_ << nl
            << "globalMeshData : sharedEdgeLabels:"
            << sharedEdgeLabelsPtr_().size() << nl
            << "globalMeshData : sharedEdgeAddr:"
            << sharedEdgeAddrPtr_().size() << endl;
    }
}


void Foam::globalMeshData::calcGlobalPointSlaves() const
{
    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointSlaves() :"
            << " calculating coupled master to slave point addressing."
            << endl;
    }

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // Calculate connected points for master points.
    globalPoints globalData(mesh_, coupledPatch(), true, true);

    globalPointSlavesPtr_.reset
    (
        new labelListList
        (
            std::move(globalData.pointPoints())
        )
    );
    globalPointTransformedSlavesPtr_.reset
    (
        new labelListList
        (
            std::move(globalData.transformedPointPoints())
        )
    );

    globalPointSlavesMapPtr_.reset
    (
        new mapDistribute
        (
            std::move(globalData.map())
        )
    );

    // Reset tag
    UPstream::msgType(oldTag);
}


void Foam::globalMeshData::calcPointConnectivity
(
    List<labelPairList>& allPointConnectivity
) const
{
    const globalIndexAndTransform& transforms = globalTransforms();
    const labelListList& slaves = globalPointSlaves();
    const labelListList& transformedSlaves = globalPointTransformedSlaves();
    const auto& slavesMap = globalPointSlavesMap();


    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // Create field with my local data
    labelPairList myData(slavesMap.constructSize());
    forAll(slaves, pointi)
    {
        myData[pointi] = transforms.encode
        (
            Pstream::myProcNo(),
            pointi,
            transforms.nullTransformIndex()
        );
    }
    // Send to master
    slavesMap.distribute(myData, true, UPstream::msgType());


    // String of connected points with their transform
    allPointConnectivity.setSize(slavesMap.constructSize());
    allPointConnectivity = labelPairList(0);

    // Pass1: do the master points since these also update local slaves
    //        (e.g. from local cyclics)
    forAll(slaves, pointi)
    {
        // Reconstruct string of connected points
        const labelList& pSlaves = slaves[pointi];
        const labelList& pTransformSlaves = transformedSlaves[pointi];

        if (pSlaves.size()+pTransformSlaves.size())
        {
            labelPairList& pConnectivity = allPointConnectivity[pointi];

            pConnectivity.setSize(1+pSlaves.size()+pTransformSlaves.size());
            label connI = 0;

            // Add myself
            pConnectivity[connI++] = myData[pointi];
            // Add untransformed points
            forAll(pSlaves, i)
            {
                pConnectivity[connI++] = myData[pSlaves[i]];
            }
            // Add transformed points.
            forAll(pTransformSlaves, i)
            {
                // Get transform from index
                label transformI = slavesMap.whichTransform
                (
                    pTransformSlaves[i]
                );
                // Add transform to connectivity
                const labelPair& n = myData[pTransformSlaves[i]];
                label proci = transforms.processor(n);
                label index = transforms.index(n);
                pConnectivity[connI++] = transforms.encode
                (
                    proci,
                    index,
                    transformI
                );
            }

            // Put back in slots
            forAll(pSlaves, i)
            {
                allPointConnectivity[pSlaves[i]] = pConnectivity;
            }
            forAll(pTransformSlaves, i)
            {
                allPointConnectivity[pTransformSlaves[i]] = pConnectivity;
            }
        }
    }


    // Pass2: see if anything is still unset (should not be the case)
    forAll(slaves, pointi)
    {
        labelPairList& pConnectivity = allPointConnectivity[pointi];

        if (pConnectivity.size() == 0)
        {
            pConnectivity.setSize(1, myData[pointi]);
        }
    }


    slavesMap.reverseDistribute
    (
        slaves.size(),
        allPointConnectivity,
        true,
        UPstream::msgType()
    );

    // Reset tag
    UPstream::msgType(oldTag);
}


void Foam::globalMeshData::calcGlobalPointEdges
(
    labelListList& globalPointEdges,
    List<labelPairList>& globalPointPoints
) const
{
    const edgeList& edges = coupledPatch().edges();
    const labelListList& pointEdges = coupledPatch().pointEdges();
    const globalIndex& globalEdgeNumbers = globalEdgeNumbering();
    const labelListList& slaves = globalPointSlaves();
    const labelListList& transformedSlaves = globalPointTransformedSlaves();
    const globalIndexAndTransform& transforms = globalTransforms();


    // Create local version
    globalPointEdges.setSize(globalPointSlavesMap().constructSize());
    globalPointPoints.setSize(globalPointSlavesMap().constructSize());
    forAll(pointEdges, pointi)
    {
        const labelList& pEdges = pointEdges[pointi];
        globalPointEdges[pointi] = globalEdgeNumbers.toGlobal(pEdges);

        labelPairList& globalPPoints = globalPointPoints[pointi];
        globalPPoints.setSize(pEdges.size());
        forAll(pEdges, i)
        {
            label otherPointi = edges[pEdges[i]].otherVertex(pointi);
            globalPPoints[i] = transforms.encode
            (
                Pstream::myProcNo(),
                otherPointi,
                transforms.nullTransformIndex()
            );
        }
    }

    // Pull slave data to master. Dummy transform.

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    globalPointSlavesMap().distribute
    (
        globalPointEdges,
        true,
        UPstream::msgType()
    );
    // Make sure second send uses 'far' away tags in case of NBX deciding on
    // multi-pass spraying of messages with consecutive tags
    globalPointSlavesMap().distribute
    (
        globalPointPoints,
        true,
        UPstream::msgType()+23456    // Unique, far enough away tag
    );

    // Add all pointEdges
    forAll(slaves, pointi)
    {
        const labelList& pSlaves = slaves[pointi];
        const labelList& pTransformSlaves = transformedSlaves[pointi];

        label n = 0;
        forAll(pSlaves, i)
        {
            n += globalPointEdges[pSlaves[i]].size();
        }
        forAll(pTransformSlaves, i)
        {
            n += globalPointEdges[pTransformSlaves[i]].size();
        }

        // Add all the point edges of the slaves to those of the (master) point
        {
            labelList& globalPEdges = globalPointEdges[pointi];
            label sz = globalPEdges.size();
            globalPEdges.setSize(sz+n);
            forAll(pSlaves, i)
            {
                const labelList& otherData = globalPointEdges[pSlaves[i]];
                forAll(otherData, j)
                {
                    globalPEdges[sz++] = otherData[j];
                }
            }
            forAll(pTransformSlaves, i)
            {
                const labelList& otherData =
                    globalPointEdges[pTransformSlaves[i]];
                forAll(otherData, j)
                {
                    globalPEdges[sz++] = otherData[j];
                }
            }

            // Put back in slots
            forAll(pSlaves, i)
            {
                globalPointEdges[pSlaves[i]] = globalPEdges;
            }
            forAll(pTransformSlaves, i)
            {
                globalPointEdges[pTransformSlaves[i]] = globalPEdges;
            }
        }


        // Same for corresponding pointPoints
        {
            labelPairList& globalPPoints = globalPointPoints[pointi];
            label sz = globalPPoints.size();
            globalPPoints.setSize(sz + n);

            // Add untransformed points
            forAll(pSlaves, i)
            {
                const labelPairList& otherData = globalPointPoints[pSlaves[i]];
                forAll(otherData, j)
                {
                    globalPPoints[sz++] = otherData[j];
                }
            }
            // Add transformed points.
            forAll(pTransformSlaves, i)
            {
                // Get transform from index
                label transformI = globalPointSlavesMap().whichTransform
                (
                    pTransformSlaves[i]
                );

                const labelPairList& otherData =
                    globalPointPoints[pTransformSlaves[i]];
                forAll(otherData, j)
                {
                    // Add transform to connectivity
                    const labelPair& n = otherData[j];
                    label proci = transforms.processor(n);
                    label index = transforms.index(n);
                    globalPPoints[sz++] = transforms.encode
                    (
                        proci,
                        index,
                        transformI
                    );
                }
            }

            // Put back in slots
            forAll(pSlaves, i)
            {
                globalPointPoints[pSlaves[i]] = globalPPoints;
            }
            forAll(pTransformSlaves, i)
            {
                globalPointPoints[pTransformSlaves[i]] = globalPPoints;
            }
        }
    }
    // Push back
    globalPointSlavesMap().reverseDistribute
    (
        slaves.size(),
        globalPointEdges,
        true,
        UPstream::msgType()
    );
    // Push back
    globalPointSlavesMap().reverseDistribute
    (
        slaves.size(),
        globalPointPoints,
        true,
        UPstream::msgType()+65432   // Unique, far enough away tag
    );

    // Reset tag
    UPstream::msgType(oldTag);
}


Foam::label Foam::globalMeshData::findTransform
(
    const labelPairList& info,
    const labelPair& remotePoint,
    const label localPoint
) const
{
    const globalIndexAndTransform& transforms = globalTransforms();

    const label remoteProci = transforms.processor(remotePoint);
    const label remoteIndex = transforms.index(remotePoint);

    label remoteTransformI = -1;
    label localTransformI = -1;
    forAll(info, i)
    {
        label proci = transforms.processor(info[i]);
        label pointi = transforms.index(info[i]);
        label transformI = transforms.transformIndex(info[i]);

        if (proci == Pstream::myProcNo() && pointi == localPoint)
        {
            localTransformI = transformI;
            //Pout<< "For local :" << localPoint
            //    << " found transform:" << localTransformI
            //    << endl;
        }
        if (proci == remoteProci && pointi == remoteIndex)
        {
            remoteTransformI = transformI;
            //Pout<< "For remote:" << remotePoint
            //    << " found transform:" << remoteTransformI
            //    << " at index:" << i
            //    << endl;
        }
    }

    if (remoteTransformI == -1 || localTransformI == -1)
    {
        FatalErrorInFunction
            << "Problem. Cannot find " << remotePoint
            << " or " << localPoint  << " "
            << coupledPatch().localPoints()[localPoint]
            << " in " << info
            << endl
            << "remoteTransformI:" << remoteTransformI << endl
            << "localTransformI:" << localTransformI
            << abort(FatalError);
    }

    return transforms.subtractTransformIndex
    (
        remoteTransformI,
        localTransformI
    );
}


void Foam::globalMeshData::calcGlobalEdgeSlaves() const
{
    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalEdgeSlaves() :"
            << " calculating coupled master to slave edge addressing." << endl;
    }

    const edgeList& edges = coupledPatch().edges();
    const globalIndex& globalEdgeNumbers = globalEdgeNumbering();
    const globalIndexAndTransform& transforms = globalTransforms();


    // The whole problem with deducting edge-connectivity from
    // point-connectivity is that one of the endpoints might be
    // a local master but the other endpoint might not. So we first
    // need to make sure that all points know about connectivity and
    // the transformations.

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // 1. collect point connectivity - basically recreating globalPoints output.
    // All points will now have a string of coupled points. The transforms are
    // in respect to the master.
    List<labelPairList> allPointConnectivity;
    calcPointConnectivity(allPointConnectivity);


    // 2. Get all pointEdges and pointPoints
    // Coupled point to global coupled edges and corresponding endpoint.
    labelListList globalPointEdges;
    List<labelPairList> globalPointPoints;
    calcGlobalPointEdges(globalPointEdges, globalPointPoints);


    // 3. Now all points have
    //      - all the connected points with original transform
    //      - all the connected global edges

    // Now all we need to do is go through all the edges and check
    // both endpoints. If there is a edge between the two which is
    // produced by transforming both points in the same way it is a shared
    // edge.

    // Collect strings of connected edges.
    List<labelPairList> allEdgeConnectivity(edges.size());

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        const labelList& pEdges0 = globalPointEdges[e[0]];
        const labelPairList& pPoints0 = globalPointPoints[e[0]];
        const labelList& pEdges1 = globalPointEdges[e[1]];
        const labelPairList& pPoints1 = globalPointPoints[e[1]];

        // Most edges will be size 2
        DynamicList<labelPair> eEdges(2);
        // Append myself.
        eEdges.append
        (
            transforms.encode
            (
                Pstream::myProcNo(),
                edgeI,
                transforms.nullTransformIndex()
            )
        );

        forAll(pEdges0, i)
        {
            forAll(pEdges1, j)
            {
                if
                (
                    pEdges0[i] == pEdges1[j]
                 && pEdges0[i] != globalEdgeNumbers.toGlobal(edgeI)
                )
                {
                    // Found a shared edge. Now check if the endpoints
                    // go through the same transformation.
                    // Local: e[0]    remote:pPoints1[j]
                    // Local: e[1]    remote:pPoints0[i]


                    // Find difference in transforms to go from point on remote
                    // edge (pPoints1[j]) to this point.

                    label transform0 = findTransform
                    (
                        allPointConnectivity[e[0]],
                        pPoints1[j],
                        e[0]
                    );
                    label transform1 = findTransform
                    (
                        allPointConnectivity[e[1]],
                        pPoints0[i],
                        e[1]
                    );

                    if (transform0 == transform1)
                    {
                        label proci = globalEdgeNumbers.whichProcID(pEdges0[i]);
                        eEdges.append
                        (
                            transforms.encode
                            (
                                proci,
                                globalEdgeNumbers.toLocal(proci, pEdges0[i]),
                                transform0
                            )
                        );
                    }
                }
            }
        }

        allEdgeConnectivity[edgeI].transfer(eEdges);
        Foam::sort
        (
            allEdgeConnectivity[edgeI],
            globalIndexAndTransform::less(transforms)
        );
    }

    // We now have - in allEdgeConnectivity - a list of edges which are shared
    // between multiple processors. Filter into non-transformed and transformed
    // connections.

    globalEdgeSlavesPtr_.reset(new labelListList(edges.size()));
    labelListList& globalEdgeSlaves = globalEdgeSlavesPtr_();
    List<labelPairList> transformedEdges(edges.size());
    forAll(allEdgeConnectivity, edgeI)
    {
        const labelPairList& edgeInfo = allEdgeConnectivity[edgeI];
        if (edgeInfo.size() >= 2)
        {
            const labelPair& masterInfo = edgeInfo[0];

            // Check if master edge (= first element (since sorted)) is me.
            if
            (
                (
                    transforms.processor(masterInfo)
                 == Pstream::myProcNo()
                )
             && (transforms.index(masterInfo) == edgeI)
            )
            {
                // Sort into transformed and untransformed
                labelList& eEdges = globalEdgeSlaves[edgeI];
                eEdges.setSize(edgeInfo.size()-1);

                labelPairList& trafoEEdges = transformedEdges[edgeI];
                trafoEEdges.setSize(edgeInfo.size()-1);

                label nonTransformI = 0;
                label transformI = 0;

                for (label i = 1; i < edgeInfo.size(); i++)
                {
                    const labelPair& info = edgeInfo[i];
                    label proci = transforms.processor(info);
                    label index = transforms.index(info);
                    label transform = transforms.transformIndex
                    (
                        info
                    );

                    if (transform == transforms.nullTransformIndex())
                    {
                        eEdges[nonTransformI++] = globalEdgeNumbers.toGlobal
                        (
                            proci,
                            index
                        );
                    }
                    else
                    {
                        trafoEEdges[transformI++] = info;
                    }
                }

                eEdges.setSize(nonTransformI);
                trafoEEdges.setSize(transformI);
            }
        }
    }


    // Construct map
    globalEdgeTransformedSlavesPtr_.reset(new labelListList());

    List<Map<label>> compactMap(Pstream::nProcs());
    globalEdgeSlavesMapPtr_.reset
    (
        new mapDistribute
        (
            globalEdgeNumbers,
            globalEdgeSlaves,

            transforms,
            transformedEdges,
            globalEdgeTransformedSlavesPtr_(),

            compactMap,
            UPstream::msgType()
        )
    );

    // Reset tag
    UPstream::msgType(oldTag);

    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalEdgeSlaves() :"
            << " coupled edges:" << edges.size()
            << " additional coupled edges:"
            << globalEdgeSlavesMapPtr_().constructSize() - edges.size()
            << endl;
    }
}


void Foam::globalMeshData::calcGlobalEdgeOrientation() const
{
    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalEdgeOrientation() :"
            << " calculating edge orientation w.r.t. master edge." << endl;
    }

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    const globalIndex& globalPoints = globalPointNumbering();

    // 1. Determine master point
    labelList masterPoint;
    {
        const mapDistribute& map = globalPointSlavesMap();

        masterPoint.setSize(map.constructSize());
        masterPoint = labelMax;

        for (label pointi = 0; pointi < coupledPatch().nPoints(); pointi++)
        {
            masterPoint[pointi] = globalPoints.toGlobal(pointi);
        }
        syncData
        (
            masterPoint,
            globalPointSlaves(),
            globalPointTransformedSlaves(),
            map,
            minEqOp<label>()
        );
    }

    // Now all points should know who is master by comparing their global
    // pointID with the masterPointID. We now can use this information
    // to find the orientation of the master edge.

    {
        const mapDistribute& map = globalEdgeSlavesMap();
        const labelListList& slaves = globalEdgeSlaves();
        const labelListList& transformedSlaves = globalEdgeTransformedSlaves();

        // Distribute orientation of master edge (in masterPoint numbering)
        labelPairList masterEdgeVerts(map.constructSize());
        masterEdgeVerts = labelPair(labelMax, labelMax);

        for (label edgeI = 0; edgeI < coupledPatch().nEdges(); edgeI++)
        {
            if
            (
                (
                    slaves[edgeI].size()
                  + transformedSlaves[edgeI].size()
                )
              > 0
            )
            {
                // I am master. Fill in my masterPoint equivalent.

                const edge& e = coupledPatch().edges()[edgeI];
                masterEdgeVerts[edgeI] = labelPair
                (
                    masterPoint[e[0]],
                    masterPoint[e[1]]
                );
            }
        }
        syncData
        (
            masterEdgeVerts,
            slaves,
            transformedSlaves,
            map,
            minEqOp<labelPair>()
        );

        // Now check my edges on how they relate to the master's edgeVerts
        globalEdgeOrientationPtr_.reset
        (
            new bitSet(coupledPatch().nEdges())
        );
        bitSet& globalEdgeOrientation = globalEdgeOrientationPtr_();

        forAll(coupledPatch().edges(), edgeI)
        {
            // Test that edge is not single edge on cyclic baffle
            if (masterEdgeVerts[edgeI] != labelPair(labelMax, labelMax))
            {
                const edge& e = coupledPatch().edges()[edgeI];
                const labelPair masterE
                (
                    masterPoint[e[0]],
                    masterPoint[e[1]]
                );

                const int stat = labelPair::compare
                (
                    masterE,
                    masterEdgeVerts[edgeI]
                );
                if (stat == 0)
                {
                    FatalErrorInFunction
                        << "problem : my edge:" << e
                        << " in master points:" << masterE
                        << " v.s. masterEdgeVerts:" << masterEdgeVerts[edgeI]
                        << exit(FatalError);
                }
                else
                {
                    globalEdgeOrientation.set(edgeI, (stat == 1));
                }
            }
            else
            {
                globalEdgeOrientation.set(edgeI, true);
            }
        }
    }

    // Reset tag
    UPstream::msgType(oldTag);

    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalEdgeOrientation() :"
            << " finished calculating edge orientation."
            << endl;
    }
}


void Foam::globalMeshData::calcPointBoundaryFaces
(
    labelListList& pointBoundaryFaces
) const
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    const Map<label>& meshPointMap = coupledPatch().meshPointMap();

    // 1. Count

    labelList nPointFaces(coupledPatch().nPoints(), Zero);

    for (const polyPatch& pp : bMesh)
    {
        if (!pp.coupled())
        {
            for (const face& f : pp)
            {
                forAll(f, fp)
                {
                    const auto iter = meshPointMap.cfind(f[fp]);
                    if (iter.good())
                    {
                        nPointFaces[iter.val()]++;
                    }
                }
            }
        }
    }


    // 2. Size

    pointBoundaryFaces.setSize(coupledPatch().nPoints());
    forAll(nPointFaces, pointi)
    {
        pointBoundaryFaces[pointi].setSize(nPointFaces[pointi]);
    }
    nPointFaces = 0;


    // 3. Fill

    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];

        if (!pp.coupled())
        {
            forAll(pp, i)
            {
                const face& f = pp[i];
                forAll(f, fp)
                {
                    const auto iter = meshPointMap.cfind(f[fp]);

                    if (iter.good())
                    {
                        label bFacei =
                             pp.start() + i - mesh_.nInternalFaces();
                        pointBoundaryFaces[iter()][nPointFaces[iter()]++] =
                            bFacei;
                    }
                }
            }
        }
    }
}


void Foam::globalMeshData::calcGlobalPointBoundaryFaces() const
{
    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointBoundaryFaces() :"
            << " calculating coupled point to boundary face addressing."
            << endl;
    }

    const label myProci = UPstream::myProcNo();

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // Construct local point to (uncoupled)boundaryfaces.
    labelListList pointBoundaryFaces;
    calcPointBoundaryFaces(pointBoundaryFaces);


    // Global indices for boundary faces
    globalBoundaryFaceNumberingPtr_.reset
    (
        new globalIndex(mesh_.nBoundaryFaces())
    );
    const auto& globalIndices = *globalBoundaryFaceNumberingPtr_;


    // Convert local boundary faces to global numbering
    globalPointBoundaryFacesPtr_.reset
    (
        new labelListList(globalPointSlavesMap().constructSize())
    );
    auto& globalPointBoundaryFaces = *globalPointBoundaryFacesPtr_;

    forAll(pointBoundaryFaces, pointi)
    {
        globalPointBoundaryFaces[pointi] = globalIndices.toGlobal
        (
            myProci,
            pointBoundaryFaces[pointi]
        );
    }


    // Pull slave pointBoundaryFaces to master
    globalPointSlavesMap().distribute
    (
        globalPointBoundaryFaces,
        true    // put data on transformed points into correct slots
    );


    // Merge slave labels into master globalPointBoundaryFaces.
    // Split into untransformed and transformed values.
    const labelListList& pointSlaves = globalPointSlaves();
    const labelListList& pointTransformSlaves =
        globalPointTransformedSlaves();
    const globalIndexAndTransform& transforms = globalTransforms();


    // Any faces coming in through transformation
    List<labelPairList> transformedFaces(pointSlaves.size());


    forAll(pointSlaves, pointi)
    {
        const labelList& slaves = pointSlaves[pointi];
        const labelList& transformedSlaves = pointTransformSlaves[pointi];

        if (slaves.size() > 0)
        {
            labelList& myBFaces = globalPointBoundaryFaces[pointi];
            label sz = myBFaces.size();

            // Count
            label n = 0;
            forAll(slaves, i)
            {
                n += globalPointBoundaryFaces[slaves[i]].size();
            }
            // Fill
            myBFaces.setSize(sz+n);
            n = sz;
            forAll(slaves, i)
            {
                const labelList& slaveBFaces =
                    globalPointBoundaryFaces[slaves[i]];

                // Add all slaveBFaces. Note that need to check for
                // uniqueness only in case of cyclics.

                for (const label slave : slaveBFaces)
                {
                    if (!SubList<label>(myBFaces, sz).found(slave))
                    {
                        myBFaces[n++] = slave;
                    }
                }
            }
            myBFaces.setSize(n);
        }


        if (transformedSlaves.size() > 0)
        {
            const labelList& untrafoFaces = globalPointBoundaryFaces[pointi];

            labelPairList& myBFaces = transformedFaces[pointi];
            label sz = myBFaces.size();

            // Count
            label n = 0;
            forAll(transformedSlaves, i)
            {
                n += globalPointBoundaryFaces[transformedSlaves[i]].size();
            }
            // Fill
            myBFaces.setSize(sz+n);
            n = sz;
            forAll(transformedSlaves, i)
            {
                label transformI = globalPointSlavesMap().whichTransform
                (
                    transformedSlaves[i]
                );

                const labelList& slaveBFaces =
                    globalPointBoundaryFaces[transformedSlaves[i]];

                for (const label slave : slaveBFaces)
                {
                    // Check that same face not already present untransformed
                    if (!untrafoFaces.found(slave))
                    {
                        label proci = globalIndices.whichProcID(slave);
                        label facei = globalIndices.toLocal(proci, slave);

                        myBFaces[n++] = transforms.encode
                        (
                            proci,
                            facei,
                            transformI
                        );
                    }
                }
            }
            myBFaces.setSize(n);
        }


        if (slaves.size() + transformedSlaves.size() == 0)
        {
             globalPointBoundaryFaces[pointi].clear();
        }
    }

    // Construct a map to get the face data directly
    List<Map<label>> compactMap(Pstream::nProcs());

    globalPointTransformedBoundaryFacesPtr_.reset
    (
        new labelListList(transformedFaces.size())
    );

    globalPointBoundaryFacesMapPtr_.reset
    (
        new mapDistribute
        (
            globalIndices,
            globalPointBoundaryFaces,

            transforms,
            transformedFaces,
            globalPointTransformedBoundaryFacesPtr_(),

            compactMap
        )
    );
    globalPointBoundaryFaces.setSize(coupledPatch().nPoints());
    globalPointTransformedBoundaryFacesPtr_().setSize(coupledPatch().nPoints());

    // Reset tag
    UPstream::msgType(oldTag);

    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointBoundaryFaces() :"
            << " coupled points:" << coupledPatch().nPoints()
            << " local boundary faces:" <<  globalIndices.localSize()
            << " additional coupled faces:"
            <<  globalPointBoundaryFacesMapPtr_().constructSize()
              - globalIndices.localSize()
            << endl;
    }
}


void Foam::globalMeshData::calcGlobalPointBoundaryCells() const
{
    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointBoundaryCells() :"
            << " calculating coupled point to boundary cell addressing."
            << endl;
    }

    const label myProci = UPstream::myProcNo();

    // Create map of boundary cells and point-cell addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label bCelli = 0;
    Map<label> meshCellMap(4*coupledPatch().nPoints());
    DynamicList<label> cellMap(meshCellMap.size());

    // Create addressing for point to boundary cells (local)
    labelListList pointBoundaryCells(coupledPatch().nPoints());

    forAll(coupledPatch().meshPoints(), pointi)
    {
        label meshPointi = coupledPatch().meshPoints()[pointi];
        const labelList& pCells = mesh_.pointCells(meshPointi);

        labelList& bCells = pointBoundaryCells[pointi];
        bCells.setSize(pCells.size());

        forAll(pCells, i)
        {
            const label celli = pCells[i];
            const auto fnd = meshCellMap.cfind(celli);

            if (fnd.good())
            {
                bCells[i] = fnd();
            }
            else
            {
                meshCellMap.insert(celli, bCelli);
                cellMap.append(celli);
                bCells[i] = bCelli;
                bCelli++;
            }
        }
    }


    boundaryCellsPtr_.reset(new labelList(std::move(cellMap)));
    labelList& boundaryCells = boundaryCellsPtr_();


    // Convert point-cells to global (boundary)cell numbers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    globalBoundaryCellNumberingPtr_.reset
    (
        new globalIndex(boundaryCells.size())
    );
    const auto& globalIndices = *globalBoundaryCellNumberingPtr_;


    // Convert local boundary cells to global numbering
    globalPointBoundaryCellsPtr_.reset
    (
        new labelListList(globalPointSlavesMap().constructSize())
    );
    auto& globalPointBoundaryCells = *globalPointBoundaryCellsPtr_;

    forAll(pointBoundaryCells, pointi)
    {
        globalPointBoundaryCells[pointi] = globalIndices.toGlobal
        (
            myProci,
            pointBoundaryCells[pointi]
        );
    }


    // Pull slave pointBoundaryCells to master
    globalPointSlavesMap().distribute
    (
        globalPointBoundaryCells,
        true    // put data on transformed points into correct slots
    );


    // Merge slave labels into master globalPointBoundaryCells
    const labelListList& pointSlaves = globalPointSlaves();
    const labelListList& pointTransformSlaves =
        globalPointTransformedSlaves();
    const globalIndexAndTransform& transforms = globalTransforms();

    List<labelPairList> transformedCells(pointSlaves.size());


    forAll(pointSlaves, pointi)
    {
        const labelList& slaves = pointSlaves[pointi];
        const labelList& transformedSlaves = pointTransformSlaves[pointi];

        if (slaves.size() > 0)
        {
            labelList& myBCells = globalPointBoundaryCells[pointi];
            label sz = myBCells.size();

            // Count
            label n = 0;
            forAll(slaves, i)
            {
                n += globalPointBoundaryCells[slaves[i]].size();
            }
            // Fill
            myBCells.setSize(sz+n);
            n = sz;
            forAll(slaves, i)
            {
                const labelList& slaveBCells =
                    globalPointBoundaryCells[slaves[i]];

                // Add all slaveBCells. Note that need to check for
                // uniqueness only in case of cyclics.

                for (const label slave : slaveBCells)
                {
                    if (!SubList<label>(myBCells, sz).found(slave))
                    {
                        myBCells[n++] = slave;
                    }
                }
            }
            myBCells.setSize(n);
        }


        if (transformedSlaves.size() > 0)
        {
            const labelList& untrafoCells = globalPointBoundaryCells[pointi];

            labelPairList& myBCells = transformedCells[pointi];
            label sz = myBCells.size();

            // Count
            label n = 0;
            forAll(transformedSlaves, i)
            {
                n += globalPointBoundaryCells[transformedSlaves[i]].size();
            }
            // Fill
            myBCells.setSize(sz+n);
            n = sz;
            forAll(transformedSlaves, i)
            {
                label transformI = globalPointSlavesMap().whichTransform
                (
                    transformedSlaves[i]
                );

                const labelList& slaveBCells =
                    globalPointBoundaryCells[transformedSlaves[i]];

                for (const label slave : slaveBCells)
                {
                    // Check that same cell not already present untransformed
                    if (!untrafoCells.found(slave))
                    {
                        label proci = globalIndices.whichProcID(slave);
                        label celli = globalIndices.toLocal(proci, slave);
                        myBCells[n++] = transforms.encode
                        (
                            proci,
                            celli,
                            transformI
                        );
                    }
                }
            }
            myBCells.setSize(n);
        }

        if (slaves.size() + transformedSlaves.size() == 0)
        {
             globalPointBoundaryCells[pointi].clear();
        }
    }

    // Construct a map to get the cell data directly
    List<Map<label>> compactMap(Pstream::nProcs());

    globalPointTransformedBoundaryCellsPtr_.reset
    (
        new labelListList(transformedCells.size())
    );

    globalPointBoundaryCellsMapPtr_.reset
    (
        new mapDistribute
        (
            globalIndices,
            globalPointBoundaryCells,

            transforms,
            transformedCells,
            globalPointTransformedBoundaryCellsPtr_(),

            compactMap
        )
    );
    globalPointBoundaryCells.setSize(coupledPatch().nPoints());
    globalPointTransformedBoundaryCellsPtr_().setSize(coupledPatch().nPoints());

    // Reset tag
    UPstream::msgType(oldTag);

    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointBoundaryCells() :"
            << " coupled points:" << coupledPatch().nPoints()
            << " local boundary cells:" <<  globalIndices.localSize()
            << " additional coupled cells:"
            <<  globalPointBoundaryCellsMapPtr_().constructSize()
              - globalIndices.localSize()
            << endl;
    }
}


void Foam::globalMeshData::calcGlobalCoPointSlaves() const
{
    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalCoPointSlaves() :"
            << " calculating coupled master to collocated"
            << " slave point addressing." << endl;
    }

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // Calculate connected points for master points.
    globalPoints globalData(mesh_, coupledPatch(), true, false);

    globalCoPointSlavesPtr_.reset
    (
        new labelListList
        (
            std::move(globalData.pointPoints())
        )
    );
    globalCoPointSlavesMapPtr_.reset
    (
        new mapDistribute
        (
            std::move(globalData.map())
        )
    );

    UPstream::msgType(oldTag);

    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalCoPointSlaves() :"
            << " finished calculating coupled master to collocated"
            << " slave point addressing." << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalMeshData::globalMeshData(const polyMesh& mesh)
:
    mesh_(mesh),
    globalMeshPointAddr_(),
    globalMeshFaceAddr_(),
    globalMeshCellAddr_(),
    processorTopology_
    (
        processorTopology::New<processorPolyPatch>
        (
            mesh.boundaryMesh(),
            mesh_.comm()
        )
    ),
    processorPatches_(),
    processorPatchIndices_(),
    processorPatchNeighbours_(),
    nGlobalPoints_(-1),
    sharedPointLabelsPtr_(nullptr),
    sharedPointAddrPtr_(nullptr),
    sharedPointGlobalLabelsPtr_(nullptr),
    nGlobalEdges_(-1),
    sharedEdgeLabelsPtr_(nullptr),
    sharedEdgeAddrPtr_(nullptr)
{
    updateMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// A non-default destructor since we had incomplete types in the header
Foam::globalMeshData::~globalMeshData()
{}


void Foam::globalMeshData::clearOut()
{
    // Point
    nGlobalPoints_ = -1;
    sharedPointLabelsPtr_.clear();
    sharedPointAddrPtr_.clear();
    sharedPointGlobalLabelsPtr_.clear();

    // Edge
    nGlobalEdges_ = -1;
    sharedEdgeLabelsPtr_.clear();
    sharedEdgeAddrPtr_.clear();

    // Coupled patch
    coupledPatchPtr_.clear();
    coupledPatchMeshEdgesPtr_.clear();
    coupledPatchMeshEdgeMapPtr_.clear();
    globalTransformsPtr_.clear();

    // Point
    globalPointNumberingPtr_.clear();
    globalPointSlavesPtr_.clear();
    globalPointTransformedSlavesPtr_.clear();
    globalPointSlavesMapPtr_.clear();

    // Edge
    globalEdgeNumberingPtr_.clear();
    globalEdgeSlavesPtr_.clear();
    globalEdgeTransformedSlavesPtr_.clear();
    globalEdgeOrientationPtr_.clear();
    globalEdgeSlavesMapPtr_.clear();

    // Face
    globalBoundaryFaceNumberingPtr_.clear();
    globalPointBoundaryFacesPtr_.clear();
    globalPointTransformedBoundaryFacesPtr_.clear();
    globalPointBoundaryFacesMapPtr_.clear();

    // Cell
    boundaryCellsPtr_.clear();
    globalBoundaryCellNumberingPtr_.clear();
    globalPointBoundaryCellsPtr_.clear();
    globalPointTransformedBoundaryCellsPtr_.clear();
    globalPointBoundaryCellsMapPtr_.clear();

    // Other: collocated points
    globalCoPointSlavesPtr_.clear();
    globalCoPointSlavesMapPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::globalMeshData::sharedPointGlobalLabels() const
{
    if (!sharedPointGlobalLabelsPtr_)
    {
        sharedPointGlobalLabelsPtr_.reset
        (
            new labelList(sharedPointLabels().size())
        );
        labelList& sharedPointGlobalLabels = sharedPointGlobalLabelsPtr_();

        IOobject addrHeader
        (
            "pointProcAddressing",
            mesh_.facesInstance()/polyMesh::meshSubDir,
            mesh_,
            IOobject::MUST_READ
        );

        if (addrHeader.typeHeaderOk<labelIOList>(true))
        {
            // There is a pointProcAddressing file so use it to get labels
            // on the original mesh
            Pout<< "globalMeshData::sharedPointGlobalLabels : "
                << "Reading pointProcAddressing" << endl;

            labelIOList pointProcAddressing(addrHeader);

            const labelList& pointLabels = sharedPointLabels();

            forAll(pointLabels, i)
            {
                // Get my mesh point
                label pointi = pointLabels[i];

                // Map to mesh point of original mesh
                sharedPointGlobalLabels[i] = pointProcAddressing[pointi];
            }
        }
        else
        {
            Pout<< "globalMeshData::sharedPointGlobalLabels :"
                << " Setting pointProcAddressing to -1" << endl;

            sharedPointGlobalLabels = -1;
        }
    }

    return *sharedPointGlobalLabelsPtr_;
}


Foam::pointField Foam::globalMeshData::sharedPoints() const
{
    // Get all processors to send their shared points to master.
    // (not very efficient)

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    pointField sharedPoints(nGlobalPoints());
    const labelList& pointAddr = sharedPointAddr();
    const labelList& pointLabels = sharedPointLabels();

    if (UPstream::master())
    {
        // Master:
        // insert my own data first
        forAll(pointLabels, i)
        {
            label sharedPointi = pointAddr[i];

            sharedPoints[sharedPointi] = mesh_.points()[pointLabels[i]];
        }

        // Receive data and insert
        for (const int proci : UPstream::subProcs())
        {
            IPstream fromProc(UPstream::commsTypes::scheduled, proci);

            labelList nbrSharedPointAddr;
            pointField nbrSharedPoints;
            fromProc >> nbrSharedPointAddr >> nbrSharedPoints;

            forAll(nbrSharedPointAddr, i)
            {
                label sharedPointi = nbrSharedPointAddr[i];

                sharedPoints[sharedPointi] = nbrSharedPoints[i];
            }
        }
    }
    else
    {
        if (UPstream::parRun())
        {
            // Send address and points
            OPstream toMaster
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo()
            );
            toMaster
                << pointAddr
                << pointField(mesh_.points(), pointLabels);
        }
    }

    // Broadcast: sharedPoints to all
    Pstream::broadcast(sharedPoints);  // == worldComm

    // Reset tag
    UPstream::msgType(oldTag);

    return sharedPoints;
}


Foam::pointField Foam::globalMeshData::geometricSharedPoints() const
{
    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // Get coords of my shared points
    pointField sharedPoints(mesh_.points(), sharedPointLabels());

    // Append from all processors, globally consistent
    Pstream::combineReduce(sharedPoints, ListOps::appendEqOp<point>());

    // Merge tolerance
    scalar tolDim = matchTol_ * mesh_.bounds().mag();

    labelList pointMap;
    Foam::inplaceMergePoints
    (
        sharedPoints,   // coordinates to merge
        tolDim,         // tolerance
        false,          // verbosity
        pointMap
    );

    // Reset tag
    UPstream::msgType(oldTag);

    return sharedPoints;
}


Foam::label Foam::globalMeshData::nGlobalPoints() const
{
    if (nGlobalPoints_ == -1)
    {
        calcSharedPoints();
    }
    return nGlobalPoints_;
}


const Foam::labelList& Foam::globalMeshData::sharedPointLabels() const
{
    if (!sharedPointLabelsPtr_)
    {
        calcSharedPoints();
    }
    return *sharedPointLabelsPtr_;
}


const Foam::labelList& Foam::globalMeshData::sharedPointAddr() const
{
    if (!sharedPointAddrPtr_)
    {
        calcSharedPoints();
    }
    return *sharedPointAddrPtr_;
}


Foam::label Foam::globalMeshData::nGlobalEdges() const
{
    if (nGlobalEdges_ == -1)
    {
        calcSharedEdges();
    }
    return nGlobalEdges_;
}


const Foam::labelList& Foam::globalMeshData::sharedEdgeLabels() const
{
    if (!sharedEdgeLabelsPtr_)
    {
        calcSharedEdges();
    }
    return *sharedEdgeLabelsPtr_;
}


const Foam::labelList& Foam::globalMeshData::sharedEdgeAddr() const
{
    if (!sharedEdgeAddrPtr_)
    {
        calcSharedEdges();
    }
    return *sharedEdgeAddrPtr_;
}


const Foam::indirectPrimitivePatch& Foam::globalMeshData::coupledPatch() const
{
    if (!coupledPatchPtr_)
    {
        const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

        label nCoupled = 0;

        forAll(bMesh, patchi)
        {
            const polyPatch& pp = bMesh[patchi];

            if (pp.coupled())
            {
                nCoupled += pp.size();
            }
        }
        labelList coupledFaces(nCoupled);
        nCoupled = 0;

        forAll(bMesh, patchi)
        {
            const polyPatch& pp = bMesh[patchi];

            if (pp.coupled())
            {
                label facei = pp.start();

                forAll(pp, i)
                {
                    coupledFaces[nCoupled++] = facei++;
                }
            }
        }

        coupledPatchPtr_.reset
        (
            new indirectPrimitivePatch
            (
                IndirectList<face>
                (
                    mesh_.faces(),
                    coupledFaces
                ),
                mesh_.points()
            )
        );

        if (debug)
        {
            Pout<< "globalMeshData::coupledPatch() :"
                << " constructed  coupled faces patch:"
                << "  faces:" << coupledPatchPtr_().size()
                << "  points:" << coupledPatchPtr_().nPoints()
                << endl;
        }
    }
    return *coupledPatchPtr_;
}


const Foam::labelList& Foam::globalMeshData::coupledPatchMeshEdges() const
{
    if (!coupledPatchMeshEdgesPtr_)
    {
        coupledPatchMeshEdgesPtr_.reset
        (
            new labelList
            (
                coupledPatch().meshEdges
                (
                    mesh_.edges(),
                    mesh_.pointEdges()
                )
            )
        );
    }
    return *coupledPatchMeshEdgesPtr_;
}


const Foam::Map<Foam::label>& Foam::globalMeshData::coupledPatchMeshEdgeMap()
const
{
    if (!coupledPatchMeshEdgeMapPtr_)
    {
        const labelList& me = coupledPatchMeshEdges();

        coupledPatchMeshEdgeMapPtr_.reset
        (
            new Map<label>(invertToMap(me))
        );
    }
    return *coupledPatchMeshEdgeMapPtr_;
}


const Foam::globalIndex& Foam::globalMeshData::globalPointNumbering() const
{
    if (!globalPointNumberingPtr_)
    {
        // Allocate unique tag for all comms
        const int oldTag = UPstream::incrMsgType();

        globalPointNumberingPtr_.reset
        (
            new globalIndex(coupledPatch().nPoints())
        );

        // Reset tag
        UPstream::msgType(oldTag);
    }
    return *globalPointNumberingPtr_;
}


const Foam::globalIndexAndTransform&
Foam::globalMeshData::globalTransforms() const
{
    if (!globalTransformsPtr_)
    {
        // Allocate unique tag for all comms
        const int oldTag = UPstream::incrMsgType();

        globalTransformsPtr_.reset(new globalIndexAndTransform(mesh_));

        // Reset tag
        UPstream::msgType(oldTag);
    }
    return *globalTransformsPtr_;
}


const Foam::labelListList& Foam::globalMeshData::globalPointSlaves() const
{
    if (!globalPointSlavesPtr_)
    {
        calcGlobalPointSlaves();
    }
    return *globalPointSlavesPtr_;
}


const Foam::labelListList& Foam::globalMeshData::globalPointTransformedSlaves()
const
{
    if (!globalPointTransformedSlavesPtr_)
    {
        calcGlobalPointSlaves();
    }
    return *globalPointTransformedSlavesPtr_;
}


const Foam::mapDistribute& Foam::globalMeshData::globalPointSlavesMap() const
{
    if (!globalPointSlavesMapPtr_)
    {
        calcGlobalPointSlaves();
    }
    return *globalPointSlavesMapPtr_;
}


const Foam::globalIndex& Foam::globalMeshData::globalEdgeNumbering() const
{
    if (!globalEdgeNumberingPtr_)
    {
        // Allocate unique tag for all comms
        const int oldTag = UPstream::incrMsgType();

        globalEdgeNumberingPtr_.reset
        (
            new globalIndex(coupledPatch().nEdges())
        );

        // Reset tag
        UPstream::msgType(oldTag);
    }
    return *globalEdgeNumberingPtr_;
}


const Foam::labelListList& Foam::globalMeshData::globalEdgeSlaves() const
{
    if (!globalEdgeSlavesPtr_)
    {
        calcGlobalEdgeSlaves();
    }
    return *globalEdgeSlavesPtr_;
}


const Foam::labelListList& Foam::globalMeshData::globalEdgeTransformedSlaves()
const
{
    if (!globalEdgeTransformedSlavesPtr_)
    {
        calcGlobalEdgeSlaves();
    }
    return *globalEdgeTransformedSlavesPtr_;
}


const Foam::bitSet& Foam::globalMeshData::globalEdgeOrientation() const
{
    if (!globalEdgeOrientationPtr_)
    {
        calcGlobalEdgeOrientation();
    }
    return *globalEdgeOrientationPtr_;
}


const Foam::mapDistribute& Foam::globalMeshData::globalEdgeSlavesMap() const
{
    if (!globalEdgeSlavesMapPtr_)
    {
        calcGlobalEdgeSlaves();
    }
    return *globalEdgeSlavesMapPtr_;
}


const Foam::globalIndex& Foam::globalMeshData::globalBoundaryFaceNumbering()
const
{
    if (!globalBoundaryFaceNumberingPtr_)
    {
        calcGlobalPointBoundaryFaces();
    }
    return *globalBoundaryFaceNumberingPtr_;
}


const Foam::labelListList& Foam::globalMeshData::globalPointBoundaryFaces()
const
{
    if (!globalPointBoundaryFacesPtr_)
    {
        calcGlobalPointBoundaryFaces();
    }
    return *globalPointBoundaryFacesPtr_;
}


const Foam::labelListList&
Foam::globalMeshData::globalPointTransformedBoundaryFaces() const
{
    if (!globalPointTransformedBoundaryFacesPtr_)
    {
        calcGlobalPointBoundaryFaces();
    }
    return *globalPointTransformedBoundaryFacesPtr_;
}


const Foam::mapDistribute& Foam::globalMeshData::globalPointBoundaryFacesMap()
const
{
    if (!globalPointBoundaryFacesMapPtr_)
    {
        calcGlobalPointBoundaryFaces();
    }
    return *globalPointBoundaryFacesMapPtr_;
}


const Foam::labelList& Foam::globalMeshData::boundaryCells() const
{
    if (!boundaryCellsPtr_)
    {
        calcGlobalPointBoundaryCells();
    }
    return *boundaryCellsPtr_;
}


const Foam::globalIndex& Foam::globalMeshData::globalBoundaryCellNumbering()
const
{
    if (!globalBoundaryCellNumberingPtr_)
    {
        calcGlobalPointBoundaryCells();
    }
    return *globalBoundaryCellNumberingPtr_;
}


const Foam::labelListList& Foam::globalMeshData::globalPointBoundaryCells()
const
{
    if (!globalPointBoundaryCellsPtr_)
    {
        calcGlobalPointBoundaryCells();
    }
    return *globalPointBoundaryCellsPtr_;
}


const Foam::labelListList&
Foam::globalMeshData::globalPointTransformedBoundaryCells() const
{
    if (!globalPointTransformedBoundaryCellsPtr_)
    {
        calcGlobalPointBoundaryCells();
    }
    return *globalPointTransformedBoundaryCellsPtr_;
}


const Foam::mapDistribute& Foam::globalMeshData::globalPointBoundaryCellsMap()
const
{
    if (!globalPointBoundaryCellsMapPtr_)
    {
        calcGlobalPointBoundaryCells();
    }
    return *globalPointBoundaryCellsMapPtr_;
}


const Foam::labelListList& Foam::globalMeshData::globalCoPointSlaves() const
{
    if (!globalCoPointSlavesPtr_)
    {
        calcGlobalCoPointSlaves();
    }
    return *globalCoPointSlavesPtr_;
}


const Foam::mapDistribute& Foam::globalMeshData::globalCoPointSlavesMap() const
{
    if (!globalCoPointSlavesMapPtr_)
    {
        calcGlobalCoPointSlaves();
    }
    return *globalCoPointSlavesMapPtr_;
}


Foam::autoPtr<Foam::globalIndex> Foam::globalMeshData::mergePoints
(
    labelList& pointToGlobal,
    labelList& uniquePoints
) const
{
    const indirectPrimitivePatch& cpp = coupledPatch();
    const globalIndex& globalCoupledPoints = globalPointNumbering();
    // Use collocated only
    const labelListList& pointSlaves = globalCoPointSlaves();
    const mapDistribute& pointSlavesMap = globalCoPointSlavesMap();

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // Points are either
    // - master with slaves
    // - slave with a master
    // - other (since e.g. non-collocated cyclics not connected)

    labelList masterGlobalPoint(cpp.nPoints(), -1);
    forAll(masterGlobalPoint, pointi)
    {
        const labelList& slavePoints = pointSlaves[pointi];
        if (slavePoints.size() > 0)
        {
            masterGlobalPoint[pointi] = globalCoupledPoints.toGlobal(pointi);
        }
    }

    // Sync by taking max
    syncData
    (
        masterGlobalPoint,
        pointSlaves,
        labelListList(0),   // no transforms
        pointSlavesMap,
        maxEqOp<label>()
    );


    // 1. Count number of masters on my processor.
    label nMaster = 0;
    bitSet isMaster(mesh_.nPoints(), true);
    forAll(pointSlaves, pointi)
    {
        if (masterGlobalPoint[pointi] == -1)
        {
            // unconnected point (e.g. from separated cyclic)
            nMaster++;
        }
        else if
        (
            masterGlobalPoint[pointi]
         == globalCoupledPoints.toGlobal(pointi)
        )
        {
            // connected master
            nMaster++;
        }
        else
        {
            // connected slave point
            isMaster.unset(cpp.meshPoints()[pointi]);
        }
    }

    label myUniquePoints = mesh_.nPoints() - cpp.nPoints() + nMaster;

    //Pout<< "Points :" << nl
    //    << "    mesh             : " << mesh_.nPoints() << nl
    //    << "    of which coupled : " << cpp.nPoints() << nl
    //    << "    of which master  : " << nMaster << nl
    //    << endl;


    // 2. Create global indexing for unique points.
    autoPtr<globalIndex> globalPointsPtr(new globalIndex(myUniquePoints));


    // 3. Assign global point numbers. Keep slaves unset.
    pointToGlobal.setSize(mesh_.nPoints());
    pointToGlobal = -1;
    uniquePoints.setSize(myUniquePoints);
    nMaster = 0;

    forAll(isMaster, meshPointi)
    {
        if (isMaster[meshPointi])
        {
            pointToGlobal[meshPointi] = globalPointsPtr().toGlobal(nMaster);
            uniquePoints[nMaster] = meshPointi;
            nMaster++;
        }
    }


    // 4. Push global index for coupled points to slaves.
    {
        labelList masterToGlobal(pointSlavesMap.constructSize(), -1);

        forAll(pointSlaves, pointi)
        {
            const labelList& slaves = pointSlaves[pointi];

            if (slaves.size() > 0)
            {
                // Duplicate master globalpoint into slave slots
                label meshPointi = cpp.meshPoints()[pointi];
                masterToGlobal[pointi] = pointToGlobal[meshPointi];
                forAll(slaves, i)
                {
                    masterToGlobal[slaves[i]] = masterToGlobal[pointi];
                }
            }
        }

        // Send back
        pointSlavesMap.reverseDistribute(cpp.nPoints(), masterToGlobal);

        // On slave copy master index into overal map.
        forAll(pointSlaves, pointi)
        {
            label meshPointi = cpp.meshPoints()[pointi];

            if (!isMaster[meshPointi])
            {
                pointToGlobal[meshPointi] = masterToGlobal[pointi];
            }
        }
    }

    // Restore tag
    UPstream::msgType(oldTag);

    return globalPointsPtr;
}


Foam::autoPtr<Foam::globalIndex> Foam::globalMeshData::mergePoints
(
    const labelUList& meshPoints,
    const Map<label>& /* unused: meshPointMap */,
    labelList& pointToGlobal,
    labelList& uniqueMeshPoints
) const
{
    const indirectPrimitivePatch& cpp = coupledPatch();
    const labelListList& pointSlaves = globalCoPointSlaves();
    const mapDistribute& pointSlavesMap = globalCoPointSlavesMap();


    // The patch points come in two variants:
    // - not on a coupled patch so guaranteed unique
    // - on a coupled patch
    // If the point is on a coupled patch the problem is that the
    // master-slave structure (globalPointSlaves etc.) assigns one of the
    // coupled points to be the master but this master point is not
    // necessarily on the patch itself! (it might just be connected to the
    // patch point via coupled patches).


    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    // Determine mapping:
    // - from patch point to coupled point (or -1)
    // - from coupled point to global patch point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const globalIndex globalPPoints(meshPoints.size());

    labelList patchToCoupled(meshPoints.size(), -1);
    labelList coupledToGlobalPatch(pointSlavesMap.constructSize(), -1);
    //label nCoupled = 0;

    // Note: loop over patch since usually smaller
    forAll(meshPoints, patchPointi)
    {
        label meshPointi = meshPoints[patchPointi];

        const auto iter = cpp.meshPointMap().cfind(meshPointi);

        if (iter.good())
        {
            patchToCoupled[patchPointi] = iter();
            coupledToGlobalPatch[iter()] = globalPPoints.toGlobal(patchPointi);
            //++nCoupled;
        }
    }

    //Pout<< "Patch:" << nl
    //    << "    points:" << meshPoints.size() << nl
    //    << "    of which on coupled patch:" << nCoupled << endl;


    // Determine master of connected points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Problem is that the coupled master might not be on the patch. So
    // work out the best patch-point master from all connected points.
    // - if the coupled master is on the patch it becomes the patch-point master
    // - else the slave with the lowest numbered patch point label

    // Get all data on master
    pointSlavesMap.distribute(coupledToGlobalPatch);
    forAll(pointSlaves, coupledPointi)
    {
        const labelList& slaves = pointSlaves[coupledPointi];

        if (slaves.size() > 0)
        {
            // I am master. What is the best candidate for patch-point master
            label masterI = labelMax;
            if (coupledToGlobalPatch[coupledPointi] != -1)
            {
                // I am master and on the coupled patch. Use me.
                masterI = coupledToGlobalPatch[coupledPointi];
            }
            else
            {
                // Get min of slaves as master.
                forAll(slaves, i)
                {
                    label slavePp = coupledToGlobalPatch[slaves[i]];
                    if (slavePp != -1 && slavePp < masterI)
                    {
                        masterI = slavePp;
                    }
                }
            }

            if (masterI != labelMax)
            {
                // Push back
                coupledToGlobalPatch[coupledPointi] = masterI;
                forAll(slaves, i)
                {
                    coupledToGlobalPatch[slaves[i]] = masterI;
                }
            }
        }
    }
    pointSlavesMap.reverseDistribute(cpp.nPoints(), coupledToGlobalPatch);


    // Generate compact numbering for master points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Now coupledToGlobalPatch is the globalIndex of the master point.
    // Now every processor can check whether they hold it and generate a
    // compact numbering.

    label nMasters = 0;
    forAll(meshPoints, patchPointi)
    {
        if (patchToCoupled[patchPointi] == -1)
        {
            nMasters++;
        }
        else
        {
            label coupledPointi = patchToCoupled[patchPointi];
            if
            (
                globalPPoints.toGlobal(patchPointi)
             == coupledToGlobalPatch[coupledPointi]
            )
            {
                // I am the master
                nMasters++;
            }
        }
    }

    autoPtr<globalIndex> globalPointsPtr(new globalIndex(nMasters));

    //Pout<< "Patch:" << nl
    //    << "    points:" << meshPoints.size() << nl
    //    << "    of which on coupled patch:" << nCoupled << nl
    //    << "    of which master:" << nMasters << endl;



    // Push back compact numbering for master point onto slaves
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    pointToGlobal.setSize(meshPoints.size());
    pointToGlobal = -1;
    uniqueMeshPoints.setSize(nMasters);

    // Sync master in global point numbering so all know the master point.
    // Initialise globalMaster to be -1 except at a globalMaster.
    labelList globalMaster(cpp.nPoints(), -1);

    nMasters = 0;
    forAll(meshPoints, patchPointi)
    {
        if (patchToCoupled[patchPointi] == -1)
        {
            uniqueMeshPoints[nMasters++] = meshPoints[patchPointi];
        }
        else
        {
            label coupledPointi = patchToCoupled[patchPointi];
            if
            (
                globalPPoints.toGlobal(patchPointi)
             == coupledToGlobalPatch[coupledPointi]
            )
            {
                globalMaster[coupledPointi] =
                    globalPointsPtr().toGlobal(nMasters);
                uniqueMeshPoints[nMasters++] = meshPoints[patchPointi];
            }
        }
    }


    // Sync by taking max
    syncData
    (
        globalMaster,
        pointSlaves,
        labelListList(0),   // no transforms
        pointSlavesMap,
        maxEqOp<label>()
    );


    // Now everyone has the master point in globalPointsPtr numbering. Fill
    // in the pointToGlobal map.
    nMasters = 0;
    forAll(meshPoints, patchPointi)
    {
        if (patchToCoupled[patchPointi] == -1)
        {
            pointToGlobal[patchPointi] = globalPointsPtr().toGlobal(nMasters++);
        }
        else
        {
            label coupledPointi = patchToCoupled[patchPointi];
            pointToGlobal[patchPointi] = globalMaster[coupledPointi];

            if
            (
                globalPPoints.toGlobal(patchPointi)
             == coupledToGlobalPatch[coupledPointi]
            )
            {
                nMasters++;
            }
        }
    }

    // Restore tag
    UPstream::msgType(oldTag);

    return globalPointsPtr;
}


void Foam::globalMeshData::movePoints(const pointField& newPoints)
{
    // Topology does not change and we don't store any geometry so nothing
    // needs to be done.
    // Only global transformations might change but this is not really
    // supported.
}


void Foam::globalMeshData::updateMesh()
{
    // Clear out old data
    clearOut();

    // Do processor patch addressing
    initProcAddr();

    scalar tolDim = matchTol_ * mesh_.bounds().mag();

    if (debug)
    {
        Pout<< "globalMeshData : merge dist:" << tolDim << endl;
    }


    const label comm = mesh_.comm();
    const label oldWarnComm = UPstream::commWarn(comm);

    if (UPstream::is_parallel(comm))
    {
        const label myProci = UPstream::myProcNo(comm);
        const label numProc = UPstream::nProcs(comm);

        // Gather all three sizes together
        labelList allSizes(3*numProc);
        {
            label* tup = allSizes.begin(3*myProci);
            tup[0] = mesh_.nPoints();
            tup[1] = mesh_.nFaces();
            tup[2] = mesh_.nCells();
        }

        UPstream::mpiAllGather(allSizes.data(), 3, comm);

        // Extract counts per mesh entity
        // TBD: check for label overflow?

        labelList counts(numProc);
        for (label proci = 0, idx = 0; proci < numProc; ++proci, idx += 3)
        {
            counts[proci] = allSizes[idx];
        }
        globalMeshPointAddr_.reset(counts);

        for (label proci = 0, idx = 1; proci < numProc; ++proci, idx += 3)
        {
            counts[proci] = allSizes[idx];
        }
        globalMeshFaceAddr_.reset(counts);

        for (label proci = 0, idx = 2; proci < numProc; ++proci, idx += 3)
        {
            counts[proci] = allSizes[idx];
        }
        globalMeshCellAddr_.reset(counts);
    }
    else
    {
        globalMeshPointAddr_.reset(globalIndex::gatherNone{}, mesh_.nPoints());
        globalMeshFaceAddr_.reset(globalIndex::gatherNone{}, mesh_.nFaces());
        globalMeshCellAddr_.reset(globalIndex::gatherNone{}, mesh_.nCells());
    }

    // Restore communicator settings
    UPstream::commWarn(oldWarnComm);

    if (debug)
    {
        Info<< "globalMeshData : Total points/faces/cells : ("
            << nTotalPoints() << ' '
            << nTotalFaces() << ' '
            << nTotalCells() << ')' << endl;
    }
}


// ************************************************************************* //
