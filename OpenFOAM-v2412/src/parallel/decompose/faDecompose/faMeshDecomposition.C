/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "faMeshDecomposition.H"
#include "Time.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "processorFaPatch.H"
#include "faMesh.H"
#include "OSspecific.H"
#include "Map.H"
#include "SLList.H"
#include "ListOps.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMeshDecomposition::distributeFaces()
{
    const word& polyMeshRegionName = mesh().name();

    Info<< "\nCalculating distribution of finiteArea faces" << endl;

    cpuTime decompositionTime;

    for (label proci = 0; proci < nProcs(); ++proci)
    {
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(proci)),
            false,  // No function objects
            false   // No extra controlDict libs
        );

        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb
            )
        );

        IOobject ioAddr
        (
            "procAddressing",
            "constant",
            polyMesh::meshSubDir,
            procFvMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );


        // faceProcAddressing (polyMesh)
        ioAddr.rename("faceProcAddressing");
        labelIOList fvFaceProcAddressing(ioAddr);

        labelHashSet faceProcAddressingHash(2*fvFaceProcAddressing.size());

        // If faMesh's fvPatch is a part of the global face zones, faces of that
        // patch will be present on all processors. Because of that, looping
        // through faceProcAddressing will decompose global faMesh faces to the
        // very last processor regardless of where fvPatch is really decomposed.
        // Since global faces which do not belong to specific processor are
        // located at the end of the faceProcAddressing, cutting it at
        // i = owner.size() will correctly decompose faMesh faces.
        // Vanja Skuric, 2016-04-21
        if (hasGlobalFaceZones_)
        {
            // owner (polyMesh)
            ioAddr.rename("owner");
            const label ownerSize = labelIOList(ioAddr).size();

            for (int i = 0; i < ownerSize; ++i)
            {
                faceProcAddressingHash.insert(fvFaceProcAddressing[i]);
            }
        }
        else
        {
            faceProcAddressingHash.insert
            (
                static_cast<labelList&>(fvFaceProcAddressing)
            );
        }

        forAll(faceLabels(), facei)
        {
            // With +1 for lookup in faceMap with flip encoding
            const label index = (faceLabels()[facei] + 1);

            if (faceProcAddressingHash.found(index))
            {
                faceToProc_[facei] = proci;
            }
        }
    }

    Info<< "\nFinished decomposition in "
        << decompositionTime.elapsedCpuTime()
        << " s" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshDecomposition::faMeshDecomposition
(
    const polyMesh& mesh,
    const label nProcessors,
    const dictionary& params
)
:
    faMesh(mesh),
    nProcs_(nProcessors),
    distributed_(false),
    hasGlobalFaceZones_(false),
    faceToProc_(nFaces()),
    procFaceLabels_(nProcs_),
    procMeshEdgesMap_(nProcs_),
    procNInternalEdges_(nProcs_, Zero),
    procPatchEdgeLabels_(nProcs_),
    procPatchPointAddressing_(nProcs_),
    procPatchEdgeAddressing_(nProcs_),
    procEdgeAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
    procBoundaryAddressing_(nProcs_),
    procPatchSize_(nProcs_),
    procPatchStartIndex_(nProcs_),
    procNeighbourProcessors_(nProcs_),
    procProcessorPatchSize_(nProcs_),
    procProcessorPatchStartIndex_(nProcs_),
    globallySharedPoints_(),
    cyclicParallel_(false)
{
    updateParameters(params);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMeshDecomposition::updateParameters
(
    const dictionary& params
)
{
    params.readIfPresent("distributed", distributed_);
    if (params.found("globalFaceZones"))
    {
        hasGlobalFaceZones_ = true;
    }
}


void Foam::faMeshDecomposition::decomposeMesh()
{
    // Decide which cell goes to which processor
    distributeFaces();

    const word& polyMeshRegionName = mesh().name();

    Info<< "\nDistributing faces to processors" << endl;

    labelList nLocalFaces(nProcs_, Zero);

    // Pass 1: determine local sizes, sanity check

    forAll(faceToProc_, facei)
    {
        const label proci = faceToProc_[facei];

        if (proci < 0 || proci >= nProcs_)
        {
            FatalErrorInFunction
                << "Invalid processor label " << proci
                << " for face " << facei << nl
                << abort(FatalError);
        }
        else
        {
            ++nLocalFaces[proci];
        }
    }

    // Adjust lengths
    forAll(nLocalFaces, proci)
    {
        procFaceAddressing_[proci].resize(nLocalFaces[proci]);
        nLocalFaces[proci] = 0;  // restart list
    }

    // Pass 2: fill in local lists
    forAll(faceToProc_, facei)
    {
        const label proci = faceToProc_[facei];
        const label localFacei = nLocalFaces[proci];
        ++nLocalFaces[proci];

        procFaceAddressing_[proci][localFacei] = facei;
    }


    // Find processor mesh faceLabels and ...

    for (label procI = 0; procI < nProcs(); procI++)
    {
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(procI)),
            false,  // No function objects
            false   // No extra controlDict libs
        );

        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb
            )
        );

        IOobject ioAddr
        (
            "procAddressing",
            "constant",
            polyMesh::meshSubDir,
            procFvMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );


        // pointProcAddressing (polyMesh)
        ioAddr.rename("pointProcAddressing");
        labelIOList fvPointProcAddressing(ioAddr);

        Map<label> fvFaceProcAddressingHash;

        {
            // faceProcAddressing (polyMesh)
            ioAddr.rename("faceProcAddressing");
            labelIOList fvFaceProcAddressing(ioAddr);

            fvFaceProcAddressingHash = invertToMap(fvFaceProcAddressing);
        }


        const labelList& curProcFaceAddressing = procFaceAddressing_[procI];

        labelList& curFaceLabels = procFaceLabels_[procI];

        curFaceLabels = labelList(curProcFaceAddressing.size(), -1);

        forAll(curProcFaceAddressing, faceI)
        {
            curFaceLabels[faceI] =
                fvFaceProcAddressingHash.find
                (
                    faceLabels()[curProcFaceAddressing[faceI]] + 1
                ).val();
        }

        // Create processor finite area mesh
        faMesh procMesh
        (
            procFvMesh,
            labelList(procFaceLabels_[procI])
        );

        const uindirectPrimitivePatch& patch = this->patch();
        const Map<label>& map = patch.meshPointMap();

        EdgeMap<label> edgesHash;

        const label nIntEdges = patch.nInternalEdges();

        for (label edgei = 0; edgei < nIntEdges; ++edgei)
        {
            edgesHash.insert(patch.edges()[edgei], edgesHash.size());
        }

        forAll(boundary(), patchi)
        {
            // Include emptyFaPatch
            const label size = boundary()[patchi].labelList::size();

            for (label edgei=0; edgei < size; ++edgei)
            {
                edgesHash.insert
                (
                    patch.edges()[boundary()[patchi][edgei]],
                    edgesHash.size()
                );
            }
        }


        const uindirectPrimitivePatch& procPatch = procMesh.patch();
        const labelList& procMeshPoints = procPatch.meshPoints();
        const edgeList& procEdges = procPatch.edges();

        labelList& curPatchPointAddressing = procPatchPointAddressing_[procI];
        curPatchPointAddressing.resize(procMeshPoints.size(), -1);

        forAll(procMeshPoints, pointi)
        {
            curPatchPointAddressing[pointi] =
                map[fvPointProcAddressing[procMeshPoints[pointi]]];
        }

        labelList& curPatchEdgeAddressing = procPatchEdgeAddressing_[procI];
        curPatchEdgeAddressing.resize(procEdges.size(), -1);

        Map<label>& curMap = procMeshEdgesMap_[procI];
        curMap.clear();
        curMap.reserve(procEdges.size());

        forAll(procEdges, edgeI)
        {
            edge curGlobalEdge(curPatchPointAddressing, procEdges[edgeI]);
            curPatchEdgeAddressing[edgeI] = edgesHash.find(curGlobalEdge).val();
        }

        forAll(curPatchEdgeAddressing, edgeI)
        {
            curMap.insert(curPatchEdgeAddressing[edgeI], edgeI);
        }

        procNInternalEdges_[procI] = procPatch.nInternalEdges();
    }


    Info << "\nDistributing edges to processors" << endl;

    // Loop through all internal edges and decide which processor they
    // belong to. First visit all internal edges.

    // set references to the original mesh
    const faBoundaryMesh& patches = boundary();
    const edgeList& edges = this->edges();
    const labelList& owner = edgeOwner();
    const labelList& neighbour = edgeNeighbour();

    // Memory management
    {
        List<SLList<label>> procEdgeList(nProcs());

        forAll(procEdgeList, procI)
        {
            for(label i=0; i<procNInternalEdges_[procI]; i++)
            {
                procEdgeList[procI].append(procPatchEdgeAddressing_[procI][i]);
            }
        }


        // Detect inter-processor boundaries

        // Neighbour processor for each subdomain
        List<SLList<label>> interProcBoundaries(nProcs());

        // Edge labels belonging to each inter-processor boundary
        List<SLList<SLList<label>>> interProcBEdges(nProcs());

        List<SLList<label>> procPatchIndex(nProcs());

        forAll(neighbour, edgeI)
        {
            if (faceToProc_[owner[edgeI]] != faceToProc_[neighbour[edgeI]])
            {
                // inter - processor patch edge found. Go through the list of
                // inside boundaries for the owner processor and try to find
                // this inter-processor patch.

                bool interProcBouFound = false;

                const label ownProc = faceToProc_[owner[edgeI]];
                const label neiProc = faceToProc_[neighbour[edgeI]];

                auto curInterProcBdrsOwnIter =
                    interProcBoundaries[ownProc].cbegin();

                auto curInterProcBEdgesOwnIter =
                    interProcBEdges[ownProc].begin();

                // WARNING: Synchronous SLList iterators

                for
                (
                    ;
                    curInterProcBdrsOwnIter.good()
                 && curInterProcBEdgesOwnIter.good();
                    ++curInterProcBdrsOwnIter,
                    ++curInterProcBEdgesOwnIter
                )
                {
                    if (curInterProcBdrsOwnIter() == neiProc)
                    {
                        // the inter - processor boundary exists. Add the face
                        interProcBouFound = true;

                        bool neighbourFound = false;

                        curInterProcBEdgesOwnIter().append(edgeI);

                        auto curInterProcBdrsNeiIter =
                            interProcBoundaries[neiProc].cbegin();

                        auto curInterProcBEdgesNeiIter =
                            interProcBEdges[neiProc].begin();

                        // WARNING: Synchronous SLList iterators

                        for
                        (
                            ;
                            curInterProcBdrsNeiIter.good()
                         && curInterProcBEdgesNeiIter.good();
                            ++curInterProcBdrsNeiIter,
                            ++curInterProcBEdgesNeiIter
                        )
                        {
                            if (curInterProcBdrsNeiIter() == ownProc)
                            {
                                // boundary found. Add the face
                                neighbourFound = true;

                                curInterProcBEdgesNeiIter().append(edgeI);
                            }

                            if (neighbourFound) break;
                        }

                        if (interProcBouFound && !neighbourFound)
                        {
                            FatalErrorIn
                                ("faDomainDecomposition::decomposeMesh()")
                                << "Inconsistency in inter - "
                                << "processor boundary lists for processors "
                                << ownProc << " and " << neiProc
                                << abort(FatalError);
                        }
                    }

                    if (interProcBouFound) break;
                }

                if (!interProcBouFound)
                {
                    // inter - processor boundaries do not exist and need to
                    // be created

                    // set the new addressing information

                    // owner
                    interProcBoundaries[ownProc].append(neiProc);
                    interProcBEdges[ownProc].append(SLList<label>(edgeI));

                    // neighbour
                    interProcBoundaries[neiProc].append(ownProc);
                    interProcBEdges[neiProc]
                        .append
                        (
                            SLList<label>(edgeI)
                        );
                }
            }
        }


        // Loop through patches. For cyclic boundaries detect inter-processor
        // edges; for all other, add edges to the edge list and remember start
        // and size of all patches.

        // for all processors, set the size of start index and patch size
        // lists to the number of patches in the mesh
        forAll(procPatchSize_, procI)
        {
            procPatchSize_[procI].setSize(patches.size());
            procPatchStartIndex_[procI].setSize(patches.size());
        }

        forAll(patches, patchI)
        {
            // Reset size and start index for all processors
            forAll(procPatchSize_, procI)
            {
                procPatchSize_[procI][patchI] = 0;
                procPatchStartIndex_[procI][patchI] =
                    procEdgeList[procI].size();
            }

            const label patchStart = patches[patchI].start();

//             if (!isA<cyclicFaPatch>(patches[patchI]))
            if (true)
            {
                // Normal patch. Add edges to processor where the face
                // next to the edge lives

                const labelListList& eF = patch().edgeFaces();

                const label size = patches[patchI].labelList::size();

                labelList patchEdgeFaces(size, -1);

                for(int eI=0; eI<size; eI++)
                {
                    patchEdgeFaces[eI] = eF[patches[patchI][eI]][0];
                }

                forAll(patchEdgeFaces, edgeI)
                {
                    const label curProc = faceToProc_[patchEdgeFaces[edgeI]];

                    // add the face
                    procEdgeList[curProc].append(patchStart + edgeI);

                    // increment the number of edges for this patch
                    procPatchSize_[curProc][patchI]++;
                }
            }
            else
            {
                // Cyclic patch special treatment

                const faPatch& cPatch = patches[patchI];

                const label cycOffset = cPatch.size()/2;

                // Set reference to faceCells for both patches
                const labelList::subList firstEdgeFaces
                (
                    cPatch.edgeFaces(),
                    cycOffset
                );

                const labelList::subList secondEdgeFaces
                (
                    cPatch.edgeFaces(),
                    cycOffset,
                    cycOffset
                );

                forAll(firstEdgeFaces, edgeI)
                {
                    if
                    (
                        faceToProc_[firstEdgeFaces[edgeI]]
                     != faceToProc_[secondEdgeFaces[edgeI]]
                    )
                    {
                        // This edge becomes an inter-processor boundary edge
                        // inter - processor patch edge found. Go through
                        // the list of inside boundaries for the owner
                        // processor and try to find this inter-processor
                        // patch.

                        cyclicParallel_ = true;

                        bool interProcBouFound = false;

                        const label ownProc =
                            faceToProc_[firstEdgeFaces[edgeI]];
                        const label neiProc =
                            faceToProc_[secondEdgeFaces[edgeI]];

                        auto curInterProcBdrsOwnIter =
                            interProcBoundaries[ownProc].cbegin();

                        auto curInterProcBEdgesOwnIter =
                            interProcBEdges[ownProc].begin();

                        // WARNING: Synchronous SLList iterators

                        for
                        (
                            ;
                            curInterProcBdrsOwnIter.good()
                         && curInterProcBEdgesOwnIter.good();
                            ++curInterProcBdrsOwnIter,
                            ++curInterProcBEdgesOwnIter
                        )
                        {
                            if (curInterProcBdrsOwnIter() == neiProc)
                            {
                                // the inter - processor boundary exists.
                                // Add the face
                                interProcBouFound = true;

                                bool neighbourFound = false;

                                curInterProcBEdgesOwnIter()
                                    .append(patchStart + edgeI);

                                auto curInterProcBdrsNeiIter
                                   = interProcBoundaries[neiProc].cbegin();

                                auto curInterProcBEdgesNeiIter =
                                    interProcBEdges[neiProc].begin();

                                // WARNING: Synchronous SLList iterators

                                for
                                (
                                    ;
                                    curInterProcBdrsNeiIter.good()
                                 && curInterProcBEdgesNeiIter.good();
                                    ++curInterProcBdrsNeiIter,
                                    ++curInterProcBEdgesNeiIter
                                )
                                {
                                    if (curInterProcBdrsNeiIter() == ownProc)
                                    {
                                        // boundary found. Add the face
                                        neighbourFound = true;

                                        curInterProcBEdgesNeiIter()
                                           .append
                                            (
                                                patchStart
                                              + cycOffset
                                              + edgeI
                                            );
                                    }

                                    if (neighbourFound) break;
                                }

                                if (interProcBouFound && !neighbourFound)
                                {
                                    FatalErrorIn
                                    (
                                        "faDomainDecomposition::decomposeMesh()"
                                    )   << "Inconsistency in inter-processor "
                                        << "boundary lists for processors "
                                        << ownProc << " and " << neiProc
                                        << " in cyclic boundary matching"
                                        << abort(FatalError);
                                }
                            }

                            if (interProcBouFound) break;
                        }

                        if (!interProcBouFound)
                        {
                            // inter - processor boundaries do not exist
                            // and need to be created

                            // set the new addressing information

                            // owner
                            interProcBoundaries[ownProc].append(neiProc);
                            interProcBEdges[ownProc]
                                .append(SLList<label>(patchStart + edgeI));

                            // neighbour
                            interProcBoundaries[neiProc].append(ownProc);
                            interProcBEdges[neiProc]
                               .append
                                (
                                    SLList<label>
                                    (
                                        patchStart
                                      + cycOffset
                                      + edgeI
                                    )
                                );
                        }
                    }
                    else
                    {
                        // This cyclic edge remains on the processor
                        label ownProc = faceToProc_[firstEdgeFaces[edgeI]];

                        // add the edge
                        procEdgeList[ownProc].append(patchStart + edgeI);

                        // increment the number of edges for this patch
                        procPatchSize_[ownProc][patchI]++;

                        // Note: I cannot add the other side of the cyclic
                        // boundary here because this would violate the order.
                        // They will be added in a separate loop below
                    }
                }

                // Ordering in cyclic boundaries is important.
                // Add the other half of cyclic edges for cyclic boundaries
                // that remain on the processor
                forAll(secondEdgeFaces, edgeI)
                {
                    if
                    (
                        faceToProc_[firstEdgeFaces[edgeI]]
                     == faceToProc_[secondEdgeFaces[edgeI]]
                    )
                    {
                        // This cyclic edge remains on the processor
                        label ownProc = faceToProc_[firstEdgeFaces[edgeI]];

                        // add the second edge
                        procEdgeList[ownProc].append
                            (patchStart + cycOffset + edgeI);

                        // increment the number of edges for this patch
                        procPatchSize_[ownProc][patchI]++;
                    }
                }
            }
        }

        // Convert linked lists into normal lists
        // Add inter-processor boundaries and remember start indices
        forAll(procEdgeList, procI)
        {
            // Get internal and regular boundary processor faces
            SLList<label>& curProcEdges = procEdgeList[procI];

            // Get reference to processor edge addressing
            labelList& curProcEdgeAddressing = procEdgeAddressing_[procI];

            labelList& curProcNeighbourProcessors =
                procNeighbourProcessors_[procI];

            labelList& curProcProcessorPatchSize =
                procProcessorPatchSize_[procI];

            labelList& curProcProcessorPatchStartIndex =
                procProcessorPatchStartIndex_[procI];

            // calculate the size
            label nEdgesOnProcessor = curProcEdges.size();

            for (const auto& bedges : interProcBEdges[procI])
            {
                nEdgesOnProcessor += bedges.size();
            }

            curProcEdgeAddressing.setSize(nEdgesOnProcessor);

            // Fill in the list. Calculate turning index.
            // Turning index will be -1 only for some edges on processor
            // boundaries, i.e. the ones where the current processor ID
            // is in the face which is a edge neighbour.
            // Turning index is stored as the sign of the edge addressing list

            label nEdges = 0;

            // Add internal and boundary edges
            // Remember to increment the index by one such that the
            // turning index works properly.
            for (const label procEdgei : curProcEdges)
            {
                curProcEdgeAddressing[nEdges] = procEdgei;
//                 curProcEdgeAddressing[nEdges] = procEdgei + 1;
                ++nEdges;
            }

            // Add inter-processor boundary edges. At the beginning of each
            // patch, grab the patch start index and size

            curProcNeighbourProcessors.setSize
            (
                interProcBoundaries[procI].size()
            );

            curProcProcessorPatchSize.setSize
            (
                interProcBoundaries[procI].size()
            );

            curProcProcessorPatchStartIndex.setSize
            (
                interProcBoundaries[procI].size()
            );

            label nProcPatches = 0;

            auto curInterProcBdrsIter =
                interProcBoundaries[procI].cbegin();

            auto curInterProcBEdgesIter =
                interProcBEdges[procI].cbegin();

            for
            (
                ;
                curInterProcBdrsIter.good()
             && curInterProcBEdgesIter.good();
                ++curInterProcBdrsIter,
                ++curInterProcBEdgesIter
            )
            {
                curProcNeighbourProcessors[nProcPatches] =
                    curInterProcBdrsIter();

                // Get start index for processor patch
                curProcProcessorPatchStartIndex[nProcPatches] = nEdges;

                label& curSize =
                    curProcProcessorPatchSize[nProcPatches];

                curSize = 0;

                // add faces for this processor boundary

                for (const label edgei : *curInterProcBEdgesIter)
                {
                    // add the edges

                    // Remember to increment the index by one such that the
                    // turning index works properly.
                    if (faceToProc_[owner[edgei]] == procI)
                    {
                        curProcEdgeAddressing[nEdges] = edgei;
//                      curProcEdgeAddressing[nEdges] = edgei + 1;
                    }
                    else
                    {
                        // turning edge
                        curProcEdgeAddressing[nEdges] = edgei;
//                      curProcEdgeAddressing[nEdges] = -(edgei + 1);
                    }

                    // increment the size
                    ++curSize;

                    ++nEdges;
                }

                ++nProcPatches;
            }
        }
    }

    Info << "\nCalculating processor boundary addressing" << endl;
    // For every patch of processor boundary, find the index of the original
    // patch. Mis-alignment is caused by the fact that patches with zero size
    // are omitted. For processor patches, set index to -1.
    // At the same time, filter the procPatchSize_ and procPatchStartIndex_
    // lists to exclude zero-size patches
    forAll(procPatchSize_, procI)
    {
        // Make a local copy of old lists
        const labelList oldPatchSizes = procPatchSize_[procI];

        const labelList oldPatchStarts = procPatchStartIndex_[procI];

        labelList& curPatchSizes = procPatchSize_[procI];

        labelList& curPatchStarts = procPatchStartIndex_[procI];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[procI];

        labelList& curBoundaryAddressing = procBoundaryAddressing_[procI];

        curBoundaryAddressing.setSize
        (
            oldPatchSizes.size()
          + curProcessorPatchSizes.size()
        );

        label nPatches = 0;

        forAll(oldPatchSizes, patchI)
        {
            //- Do not suppress zero sized patches since make parallel
            //  actions inside patches near impossible.
            //if (oldPatchSizes[patchI] > 0)
            {
                curBoundaryAddressing[nPatches] = patchI;

                curPatchSizes[nPatches] = oldPatchSizes[patchI];

                curPatchStarts[nPatches] = oldPatchStarts[patchI];

                nPatches++;
            }
        }

        // reset to the size of live patches
        curPatchSizes.setSize(nPatches);
        curPatchStarts.setSize(nPatches);

        forAll(curProcessorPatchSizes, procPatchI)
        {
            curBoundaryAddressing[nPatches] = -1;

            nPatches++;
        }

        curBoundaryAddressing.setSize(nPatches);
    }


    // Gather data about globally shared points

    labelList globallySharedPoints_(0);

    // Memory management
    {
        labelList pointsUsage(nPoints(), Zero);

        // Globally shared points are the ones used by more than 2 processors
        // Size the list approximately and gather the points
        labelHashSet gSharedPoints
        (
            min(100, nPoints()/1000)
        );

        // Loop through all the processors and mark up points used by
        // processor boundaries.  When a point is used twice, it is a
        // globally shared point

        for (label procI = 0; procI < nProcs(); procI++)
        {
            // Get list of edge labels
            const labelList& curEdgeLabels = procEdgeAddressing_[procI];

            // Get start of processor faces
            const labelList& curProcessorPatchStarts =
                procProcessorPatchStartIndex_[procI];

            const labelList& curProcessorPatchSizes =
                procProcessorPatchSize_[procI];

            // Reset the lookup list
            pointsUsage = 0;

            forAll(curProcessorPatchStarts, patchI)
            {
                const label curStart = curProcessorPatchStarts[patchI];
                const label curEnd = curStart + curProcessorPatchSizes[patchI];

                for
                (
                    label edgeI = curStart;
                    edgeI < curEnd;
                    edgeI++
                )
                {
                    // Mark the original edge as used
                    // Remember to decrement the index by one (turning index)
                    const label curE = curEdgeLabels[edgeI];

                    const edge& e = edges[curE];

                    forAll(e, pointI)
                    {
                        if (pointsUsage[e[pointI]] == 0)
                        {
                            // Point not previously used
                            pointsUsage[e[pointI]] = patchI + 1;
                        }
                        else if (pointsUsage[e[pointI]] != patchI + 1)
                        {
                            // Point used by some other patch = global point!
                            gSharedPoints.insert(e[pointI]);
                        }
                    }
                }
            }
        }

        // Grab the result from the hash list
        globallySharedPoints_ = gSharedPoints.sortedToc();
    }


    // Edge label for faPatches

    for (label procI = 0; procI < nProcs(); procI++)
    {
        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(procI)),
            false,  // No function objects
            false   // No extra controlDict libs
        );


        // Read volume mesh
        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb
            )
        );

        // create finite area mesh
        faMesh procMesh
        (
            procFvMesh,
            labelList(procFaceLabels_[procI])
        );


        const labelList& curEdgeAddressing = procEdgeAddressing_[procI];

        const labelList& curPatchStartIndex = procPatchStartIndex_[procI];
        const labelList& curPatchSize = procPatchSize_[procI];

        const labelList& curProcessorPatchStartIndex =
            procProcessorPatchStartIndex_[procI];

        const labelList& curProcessorPatchSize =
            procProcessorPatchSize_[procI];

        labelListList& curPatchEdgeLabels = procPatchEdgeLabels_[procI];
        curPatchEdgeLabels.resize
        (
            curPatchSize.size()
          + curProcessorPatchSize.size()
        );

        forAll(curPatchSize, patchI)
        {
            labelList& curEdgeLabels = curPatchEdgeLabels[patchI];
            curEdgeLabels.setSize(curPatchSize[patchI], -1);

            label edgeI = 0;

            for
            (
                int i=curPatchStartIndex[patchI];
                i<(curPatchStartIndex[patchI]+curPatchSize[patchI]);
                i++
            )
            {
                curEdgeLabels[edgeI] =
                    procMeshEdgesMap_[procI][curEdgeAddressing[i]];
                edgeI++;
            }
        }

        forAll(curProcessorPatchSize, patchI)
        {
            labelList& curEdgeLabels  =
                curPatchEdgeLabels[curPatchSize.size() + patchI];
            curEdgeLabels.setSize(curProcessorPatchSize[patchI], -1);

            label edgeI = 0;

            for
            (
                int i=curProcessorPatchStartIndex[patchI];
                i<(curProcessorPatchStartIndex[patchI]
                +curProcessorPatchSize[patchI]);
                i++
            )
            {
                curEdgeLabels[edgeI] =
                    procMeshEdgesMap_[procI][curEdgeAddressing[i]];
                edgeI++;
            }
        }
    }
}


bool Foam::faMeshDecomposition::writeDecomposition()
{
    const word& polyMeshRegionName = mesh().name();

    Info<< "\nConstructing processor FA meshes" << endl;

    // Make a lookup map for globally shared points
    Map<label> sharedPointLookup(invertToMap(globallySharedPoints_));


    label maxProcFaces = 0;
    label totProcFaces = 0;
    label maxProcEdges = 0;
    label totProcEdges = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;

    // Write out the meshes
    for (label procI = 0; procI < nProcs(); procI++)
    {
        // Create processor mesh without a boundary

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(procI)),
            false,  // No function objects
            false   // No extra controlDict libs
        );

        // Read volume mesh
        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb
            )
        );

        labelIOList fvBoundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                "constant",
                polyMesh::meshSubDir,
                procFvMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );


        // Create finite area mesh
        faMesh procMesh
        (
            procFvMesh,
            labelList(procFaceLabels_[procI])
        );

        // Create processor boundary patches
        const labelList& curBoundaryAddressing =
            procBoundaryAddressing_[procI];

        const labelList& curPatchSizes = procPatchSize_[procI];

        const labelList& curNeighbourProcessors =
            procNeighbourProcessors_[procI];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[procI];

        const labelListList& curPatchEdgeLabels =
            procPatchEdgeLabels_[procI];

        const faPatchList& meshPatches = boundary();

        faPatchList procPatches
        (
            curPatchSizes.size() + curProcessorPatchSizes.size()
        );

        label nPatches = 0;

        forAll(curPatchSizes, patchi)
        {
            const labelList& curEdgeLabels = curPatchEdgeLabels[nPatches];

            const label neiPolyPatchId =
                fvBoundaryProcAddressing.find
                (
                    meshPatches[curBoundaryAddressing[patchi]]
                    .ngbPolyPatchIndex()
                );

            procPatches.set
            (
                nPatches,
                meshPatches[curBoundaryAddressing[patchi]].clone
                (
                    procMesh.boundary(),
                    curEdgeLabels,
                    nPatches,
                    neiPolyPatchId
                )
            );
            ++nPatches;
        }

        forAll(curProcessorPatchSizes, procPatchI)
        {
            const labelList& curEdgeLabels = curPatchEdgeLabels[nPatches];

            procPatches.set
            (
                nPatches,
                new processorFaPatch
                (
                    curEdgeLabels,
                    nPatches,
                    procMesh.boundary(),
                    -1,
                    procI,
                    curNeighbourProcessors[procPatchI]
                )
            );

            ++nPatches;
        }

        // Add boundary patches
        procMesh.addFaPatches(procPatches);

        // More precision (for points data)
        IOstream::minPrecision(10);

        procMesh.write();

        // Statistics
        Info<< nl << "Processor " << procI;

        if (procMesh.nFaces())
        {
            Info<< nl << "    ";
        }
        else
        {
            Info<< ": ";
        }

        Info<< "Number of faces = " << procMesh.nFaces() << nl;

        if (procMesh.nFaces())
        {
            Info<< "    Number of points = " << procMesh.nPoints() << nl;
        }

        totProcFaces += procMesh.nFaces();
        maxProcFaces = max(maxProcFaces, procMesh.nFaces());

        label nBoundaryEdges = 0;
        label nProcPatches = 0;
        label nProcEdges = 0;

        for (const faPatch& fap : procMesh.boundary())
        {
            const auto* ppp = isA<processorFaPatch>(fap);

            if (ppp)
            {
                const auto& procPatch = *ppp;

                Info<< "    Number of edges shared with processor "
                    << procPatch.neighbProcNo() << " = "
                    << procPatch.size() << nl;

                nProcEdges += procPatch.size();
                ++nProcPatches;
            }
            else
            {
                nBoundaryEdges += fap.size();
            }
        }

        if (procMesh.nFaces() && (nBoundaryEdges || nProcEdges))
        {
            Info<< "    Number of processor patches = " << nProcPatches << nl
                << "    Number of processor edges = " << nProcEdges << nl
                << "    Number of boundary edges = " << nBoundaryEdges << nl;
        }

        totProcEdges += nProcEdges;
        totProcPatches += nProcPatches;
        maxProcEdges = max(maxProcEdges, nProcEdges);
        maxProcPatches = max(maxProcPatches, nProcPatches);

        // Write the addressing information

        IOobject ioAddr
        (
            "procAddressing",
            "constant",
            faMesh::meshSubDir,
            procMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        // pointProcAddressing
        ioAddr.rename("pointProcAddressing");
        IOListRef<label>(ioAddr, procPatchPointAddressing_[procI]).write();

        // edgeProcAddressing
        ioAddr.rename("edgeProcAddressing");
        IOListRef<label>(ioAddr, procEdgeAddressing_[procI]).write();

        // faceProcAddressing
        ioAddr.rename("faceProcAddressing");
        IOListRef<label>(ioAddr, procFaceAddressing_[procI]).write();

        // boundaryProcAddressing
        ioAddr.rename("boundaryProcAddressing");
        IOListRef<label>(ioAddr, procBoundaryAddressing_[procI]).write();
    }


    // Summary stats
    Info<< nl
        << "Number of processor edges = " << (totProcEdges/2) << nl
        << "Max number of faces = " << maxProcFaces;

    if (maxProcFaces != totProcFaces)
    {
        scalar avgValue = scalar(totProcFaces)/nProcs_;

        Info<< " (" << 100.0*(maxProcFaces-avgValue)/avgValue
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

    Info<< "Max number of edges between processors = " << maxProcEdges;
    if (totProcEdges)
    {
        scalar avgValue = scalar(totProcEdges)/nProcs_;

        Info<< " (" << 100.0*(maxProcEdges-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl << endl;

    return true;
}


// ************************************************************************* //
