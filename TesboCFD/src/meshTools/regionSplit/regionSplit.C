/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "regionSplit.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "globalIndex.H"
#include "syncTools.H"
#include "clockValue.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionSplit, 0);
}

static constexpr Foam::label UNASSIGNED = -1;
static constexpr Foam::label BLOCKED = -2;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- The sizes of a List of containers (eg, labelHashSet)
template<class Container>
static labelList containerSizes(const UList<Container>& input)
{
    const label len = input.size();

    labelList output(len);

    for (label i = 0; i < len; ++i)
    {
        output[i] = input[i].size();
    }

    return output;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionSplit::checkBoundaryFaceSync
(
    const boolList& blockedFace
) const
{
    if (blockedFace.size())
    {
        // Check that blockedFace is synced.
        boolList syncBlockedFace(blockedFace);
        syncTools::swapFaceList(mesh(), syncBlockedFace);

        forAll(syncBlockedFace, facei)
        {
            if
            (
                blockedFace.test(facei)
             != syncBlockedFace.test(facei)
            )
            {
                FatalErrorInFunction
                    << "Face " << facei << " not synchronised. My value:"
                    << blockedFace.test(facei) << "  coupled value:"
                    << syncBlockedFace.test(facei) << nl
                    << abort(FatalError);
            }
        }
    }
}


void Foam::regionSplit::updateFacePair
(
    const label face0,
    const label face1,

    labelList& faceRegion,
    DynamicList<label>& facesChanged
) const
{
    if (faceRegion[face0] == UNASSIGNED)
    {
        if (faceRegion[face1] >= 0)
        {
            faceRegion[face0] = faceRegion[face1];
            facesChanged.append(face0);
        }
        else if (faceRegion[face1] == BLOCKED)
        {
            // face1 blocked but not face0.
            // - illegal for coupled faces, OK for explicit connections.
        }
    }
    else if (faceRegion[face0] >= 0)
    {
        if (faceRegion[face1] == UNASSIGNED)
        {
            faceRegion[face1] = faceRegion[face0];
            facesChanged.append(face1);
        }
        else if (faceRegion[face1] == BLOCKED)
        {
            // face1 blocked but not face0.
            // - illegal for coupled faces, OK for explicit connections.
        }
        else if (faceRegion[face1] != faceRegion[face0])
        {
            FatalErrorInFunction
                << "Problem : coupled face " << face0
                << " on patch " << mesh().boundaryMesh().whichPatch(face0)
                << " has region " << faceRegion[face0]
                << " but coupled face " << face1
                << " has region " << faceRegion[face1] << nl
                << "Is your blocked faces specification"
                << " synchronized across coupled boundaries?" << endl
                << abort(FatalError);
        }
    }
}


void Foam::regionSplit::fillSeedMask
(
    const UList<labelPair>& explicitConnections,
    const label seedCellId,
    const label markValue,
    labelList& cellRegion,
    labelList& faceRegion
) const
{
    // Seed cell
    cellRegion[seedCellId] = markValue;

    // Faces on seed cell
    changedFaces_.clear();
    for (const label facei : mesh().cells()[seedCellId])
    {
        if (faceRegion[facei] == UNASSIGNED)
        {
            faceRegion[facei] = markValue;
            changedFaces_.append(facei);
        }
    }

    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    // Loop over changed faces. FaceCellWave in small.

    while (changedFaces_.size())
    {
        changedCells_.clear();

        for (const label facei : changedFaces_)
        {
            const label own = mesh().faceOwner()[facei];

            if (cellRegion[own] == UNASSIGNED)
            {
                cellRegion[own] = markValue;
                changedCells_.append(own);
            }

            if (mesh().isInternalFace(facei))
            {
                const label nei = mesh().faceNeighbour()[facei];

                if (cellRegion[nei] == UNASSIGNED)
                {
                    cellRegion[nei] = markValue;
                    changedCells_.append(nei);
                }
            }
        }

        if (debug & 2)
        {
            Pout<< " Changed cells / faces : "
                << changedCells_.size() << " / " << changedFaces_.size()
                << " before sync" << endl;
        }

        // Loop over changedCells and collect faces
        changedFaces_.clear();
        for (const label celli : changedCells_)
        {
            for (const label facei : mesh().cells()[celli])
            {
                if (faceRegion[facei] == UNASSIGNED)
                {
                    faceRegion[facei] = markValue;
                    changedFaces_.append(facei);
                }
            }
        }


        // Update locally coupled faces
        // Global connections are done later.

        for (const polyPatch& pp : patches)
        {
            const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

            if (cpp && cpp->owner())
            {
                // Transfer from neighbourPatch to here or vice versa.
                const auto& cycPatch = *cpp;

                label face0 = cycPatch.start();

                forAll(cycPatch, i)
                {
                    const label face1 = cycPatch.transformGlobalFace(face0);

                    updateFacePair
                    (
                        face0,
                        face1,
                        faceRegion,
                        changedFaces_
                    );

                    ++face0;
                }
            }
        }

        for (const labelPair& pr : explicitConnections)
        {
            updateFacePair
            (
                pr.first(),
                pr.second(),
                faceRegion,
                changedFaces_
            );
        }

        if (debug & 2)
        {
            Pout<< " Changed faces : "
                << changedFaces_.size()
                << " after sync" << endl;
        }
    }
}


Foam::label Foam::regionSplit::localRegionSplit
(
    const UList<labelPair>& explicitConnections,

    labelList& cellRegion,
    labelList& faceRegion
) const
{
    clockValue timing(debug);

    changedCells_.reserve(mesh_.nCells());
    changedFaces_.reserve(mesh_.nFaces());


    // Assign local regions
    // ~~~~~~~~~~~~~~~~~~~~

    // Start with region 0
    label nLocalRegions = 0;

    for (label seedCellId = 0; seedCellId < cellRegion.size(); ++seedCellId)
    {
        // Find next unset cell - use as seed

        for (; seedCellId < cellRegion.size(); ++seedCellId)
        {
            if (cellRegion[seedCellId] == UNASSIGNED)
            {
                break;
            }
        }

        if (seedCellId >= cellRegion.size())
        {
            break;
        }

        fillSeedMask
        (
            explicitConnections,
            seedCellId,
            nLocalRegions,
            cellRegion,
            faceRegion
        );

        ++nLocalRegions; // Next region
    }

    // Discard temporary working data
    changedCells_.clearStorage();
    changedFaces_.clearStorage();

    if (debug)
    {
        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] < 0)
            {
                FatalErrorInFunction
                    << "cell:" << celli << " region:" << cellRegion[celli]
                    << abort(FatalError);
            }
        }

        forAll(faceRegion, facei)
        {
            if (faceRegion[facei] == UNASSIGNED)
            {
                FatalErrorInFunction
                    << "face:" << facei << " region:" << faceRegion[facei]
                    << abort(FatalError);
            }
        }
    }

    DebugInfo << "regionSplit = " << double(timing.elapsed()) << "s\n";

    return nLocalRegions;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSplit::regionSplit
(
    const polyMesh& mesh,
    const bool doGlobalRegions
)
:
    regionSplit
    (
        mesh,
        bitSet(),           // No blockedFace
        List<labelPair>(),  // No explicitConnections
        doGlobalRegions
    )
{}


Foam::regionSplit::regionSplit
(
    const polyMesh& mesh,
    const bitSet& blockedFace,
    const List<labelPair>& explicitConnections,
    const bool doGlobalRegions
)
:
    MeshObject_type(mesh),
    labelList(mesh.nCells(), UNASSIGNED),
    globalNumbering_()
{
    // if (debug)
    // {
    //     checkBoundaryFaceSync(blockedFace);
    // }

    labelList& cellRegion = *this;

    labelList faceRegion(mesh.nFaces(), UNASSIGNED);

    for (const label facei : blockedFace)
    {
        faceRegion[facei] = BLOCKED;
    }

    const label numLocalRegions =
        localRegionSplit(explicitConnections, cellRegion, faceRegion);

    faceRegion.clear();

    if (doGlobalRegions)
    {
        // Wrap bitset or bools
        bitSetOrBoolList hasBlockedFace(blockedFace);

        globalNumbering_ =
            reduceRegionsImpl(numLocalRegions, hasBlockedFace, cellRegion);

    }
    else
    {
        globalNumbering_ = globalIndex(numLocalRegions);
    }
}


Foam::regionSplit::regionSplit
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections,
    const bool doGlobalRegions
)
:
    MeshObject_type(mesh),
    labelList(mesh.nCells(), UNASSIGNED),
    globalNumbering_()
{
    if (debug)
    {
        checkBoundaryFaceSync(blockedFace);
    }

    labelList& cellRegion = *this;

    labelList faceRegion(mesh.nFaces(), UNASSIGNED);

    forAll(blockedFace, facei)
    {
        if (blockedFace.test(facei))
        {
            faceRegion[facei] = BLOCKED;
        }
    }


    const label numLocalRegions =
        localRegionSplit(explicitConnections, cellRegion, faceRegion);

    faceRegion.clear();

    if (doGlobalRegions)
    {
        // Wrap bitset or bools
        bitSetOrBoolList hasBlockedFace(blockedFace);

        globalNumbering_ =
            reduceRegionsImpl(numLocalRegions, hasBlockedFace, cellRegion);
    }
    else
    {
        globalNumbering_ = globalIndex(numLocalRegions);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::globalIndex
Foam::regionSplit::reduceRegionsImpl
(
    const label numLocalRegions,
    const bitSetOrBoolList& blockedFace,
    labelList& cellRegion
) const
{
    clockValue timing(debug);

    if (cellRegion.size() != mesh().nCells())
    {
        FatalErrorInFunction
            << "The cellRegion size " << cellRegion.size()
            << " != number of cells " << mesh().nCells() << endl
            << abort(FatalError);
    }


    // (numLocalRegions < 0) to signal that region information should be
    // determined ourselves. This is not really efficient, but can be useful

    const label nLocalRegions =
    (
        numLocalRegions < 0
      ? labelHashSet(cellRegion).size()
      : numLocalRegions
    );


    // Preliminary global region numbers
    const globalIndex globalRegions(nLocalRegions);


    // Lookup table of local region to global region.
    // Initially an identity mapping of the uncombined global values.

    Map<label> localToGlobal(2*nLocalRegions);
    for (const label regioni : cellRegion)
    {
        localToGlobal.insert(regioni, globalRegions.toGlobal(regioni));
    }

    // To update the localToGlobal mapping during traversal of the boundaries
    // and later when finalizing things.
    Map<label> updateLookup(2*nLocalRegions);


    // Note that we use two separate maps during the process.
    // The localToGlobal is used to map the local to global regions.
    // Merging across processors will normally make this a many->few mapping.
    // However, we may need to walk up and down processor boundaries several
    // times before all the information propagates through.
    // During these traversals, it will normally be more efficient to just
    // update the mapping without updating the cellRegion immediately.
    // Only after everything is finalized do we renumber all of the cell
    // regions.


    // Merge global regions
    // ~~~~~~~~~~~~~~~~~~~~
    // Regions across non-blocked proc patches get merged.
    // This will set merged global regions to be the min of both.
    // (this will create gaps in the global region list so they will get
    // merged later on)

    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    // Buffer for swapping boundary information
    labelList nbrRegion(mesh().nBoundaryFaces());

    bool emitWarning = true;

    do
    {
        if (debug)
        {
            Pout<< nl << "-- Starting Iteration --" << endl;
        }

        updateLookup.clear();
        nbrRegion = UNASSIGNED;

        // Region information to send
        for (const polyPatch& pp : patches)
        {
            if (pp.coupled())
            {
                SubList<label> patchNbrRegion
                (
                    nbrRegion,
                    pp.size(),
                    pp.offset()
                );

                const labelUList& faceCells = pp.faceCells();
                forAll(faceCells, patchFacei)
                {
                    const label celli = faceCells[patchFacei];
                    const label meshFacei = pp.start()+patchFacei;

                    if (!blockedFace.test(meshFacei))
                    {
                        // Send the most currently updated region Id
                        const label orig = cellRegion[celli];

                        patchNbrRegion[patchFacei] = localToGlobal[orig];
                    }
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh(), nbrRegion);

        // Receive and reduce region information
        for (const polyPatch& pp : patches)
        {
            if (pp.coupled())
            {
                SubList<label> patchNbrRegion
                (
                    nbrRegion,
                    pp.size(),
                    pp.offset()
                );

                const labelUList& faceCells = pp.faceCells();
                forAll(faceCells, patchFacei)
                {
                    const label celli = faceCells[patchFacei];
                    const label meshFacei = pp.start()+patchFacei;

                    if (!blockedFace.test(meshFacei))
                    {
                        // Reduction by retaining the min region id.

                        const label orig = cellRegion[celli];

                        const label sent = localToGlobal[orig];
                        const label recv = patchNbrRegion[patchFacei];

                        if (recv == UNASSIGNED)
                        {
                            if (emitWarning)
                            {
                                Pout<<"Warning in regionSplit:"
                                    " received unassigned on "
                                    << pp.name() << " at patchFace "
                                    << patchFacei
                                    << ". Check synchronisation in caller"
                                    << nl;
                            }
                        }
                        else if (recv < sent)
                        {
                            // Record the minimum value seen

                            auto fnd = updateLookup.find(sent);
                            if (!fnd.good())
                            {
                                updateLookup.insert(sent, recv);
                            }
                            else if (recv < *fnd)
                            {
                                *fnd = recv;
                            }
                        }
                    }
                }
            }
        }


        // Note: by always using the minimum region number across the
        // processor faces, we effect a consolidation of connected regions
        // and converge to a unique number for each distinct region.


        // Update localToGlobal according to newly exchanged information

        inplaceMapValue(updateLookup, localToGlobal);

        if (debug & 2)
        {
            labelList keys(localToGlobal.sortedToc());
            labelList vals(keys.size());
            forAll(keys, i)
            {
                vals[i] = localToGlobal[keys[i]];
            }

            Pout<< "Updated local regions:" << nl
                << "old: " << flatOutput(keys) << nl
                << "new: " << flatOutput(vals) << endl;
        }
        else if (debug)
        {
            Pout<< "Updated " << localToGlobal.size()
                << " local regions" << endl;
        }

        emitWarning = false;
        // Continue until there are no further changes
    }
    while (returnReduceOr(updateLookup.size()));


    //
    // We will need compact numbers locally and non-locally
    //

    // Determine the local compact numbering
    label nCompact = 0;
    {
        labelHashSet localRegion(2*localToGlobal.size());

        forAllConstIters(localToGlobal, iter)
        {
            const label regioni = iter.val();

            if (globalRegions.isLocal(regioni))
            {
                localRegion.insert(regioni);
            }
        }

        nCompact = localRegion.size();
    }


    // The new global numbering using compacted local regions
    globalIndex globalCompact(nCompact);


    // Determine the following:
    // - the local compact regions (store as updateLookup)
    // - the non-local regions, ordered according to the processor on which
    //   they are local.


    // The local compaction map (updated local to compact local numbering)
    updateLookup.clear();

    labelListList sendNonLocal(Pstream::nProcs());

    {
        List<labelHashSet> nonLocal(Pstream::nProcs());

        // Use estimate of sizing for non-local regions
        label estimatedCount((nLocalRegions-nCompact)/Pstream::nProcs());

        forAll(nonLocal, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                nonLocal[proci].reserve(estimatedCount);
            }
        }


        forAllConstIters(localToGlobal, iter)
        {
            const label regioni = iter.val();

            if (globalRegions.isLocal(regioni))
            {
                updateLookup.insert
                (
                    regioni,
                    globalCompact.toGlobal(updateLookup.size())
                );
            }
            else
            {
                nonLocal[globalRegions.whichProcID(regioni)].insert(regioni);
            }
        }

        if (debug)
        {
            Pout<< " per processor nonLocal regions: "
                << flatOutput(containerSizes(nonLocal)) << endl;
        }


        // Convert to label list
        forAll(sendNonLocal, proci)
        {
            sendNonLocal[proci] = nonLocal[proci].toc();
        }
    }


    // Get the wanted region labels into recvNonLocal
    labelListList recvNonLocal;
    Pstream::exchange<labelList, label>
    (
        sendNonLocal,
        recvNonLocal
    );


    // The recvNonLocal[proci] region labels are what proci requires.
    // Transcribe into their compacted number.

    {
        labelListList sendLocal(std::move(recvNonLocal));

        for (labelList& send : sendLocal)
        {
            for (label& regioni : send)
            {
                regioni = updateLookup[regioni];
            }
        }

        // Send back (into recvNonLocal)
        Pstream::exchange<labelList, label>
        (
            sendLocal,
            recvNonLocal
        );
    }


    // Now recvNonLocal and sendNonLocal contain matched pairs with
    // sendNonLocal being the non-compact region and recvNonLocal being
    // the compact region.
    //
    // Insert these into the local compaction map.

    forAll(recvNonLocal, proci)
    {
        const labelList& send = sendNonLocal[proci];
        const labelList& recv = recvNonLocal[proci];

        forAll(send, i)
        {
            updateLookup.insert(send[i], recv[i]);
        }
    }


    // Now renumber the localToGlobal to use the final compact global values
    inplaceMapValue(updateLookup, localToGlobal);


    // Can now finally use localToGlobal to renumber cellRegion

    forAll(cellRegion, celli)
    {
        cellRegion[celli] = localToGlobal[cellRegion[celli]];
    }

    DebugInfo
        <<"regionSplit::reduceRegions = " << double(timing.elapsed()) << "s\n";

    return globalCompact;
}


Foam::globalIndex
Foam::regionSplit::reduceRegions
(
    const label numLocalRegions,
    const bitSet& blockedFace,

    labelList& cellRegion
) const
{
    // Wrap bitset or bools
    bitSetOrBoolList hasBlockedFace(blockedFace);

    return reduceRegionsImpl(numLocalRegions, hasBlockedFace, cellRegion);
}


// ************************************************************************* //
