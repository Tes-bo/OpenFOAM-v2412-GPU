/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "FaceCellWave.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "UIPstream.H"
#include "UOPstream.H"
#include "PstreamReduceOps.H"
#include "debug.H"
#include "typeInfo.H"
#include "SubField.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<class Type, class TrackingData>
    class combine
    {
        //- Combine operator for AMIInterpolation

        FaceCellWave<Type, TrackingData>& solver_;

        const cyclicAMIPolyPatch& patch_;

    public:

            combine
            (
                FaceCellWave<Type, TrackingData>& solver,
                const cyclicAMIPolyPatch& patch
            )
            :
                solver_(solver),
                patch_(patch)
            {}


            void operator()
            (
                Type& x,
                const label facei,
                const Type& y,
                const scalar weight
            ) const
            {
                if (y.valid(solver_.data()))
                {
                    label meshFacei = -1;
                    if (patch_.owner())
                    {
                        meshFacei = patch_.start() + facei;
                    }
                    else
                    {
                        meshFacei = patch_.neighbPatch().start() + facei;
                    }
                    x.updateFace
                    (
                        solver_.mesh(),
                        meshFacei,
                        y,
                        solver_.propagationTol(),
                        solver_.data()
                    );
                }
            }
    };
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, class TrackingData>
bool Foam::FaceCellWave<Type, TrackingData>::updateCell
(
    const label celli,
    const label neighbourFacei,
    const Type& neighbourInfo,
    const scalar tol,
    Type& cellInfo
)
{
    // Update info for celli, at position pt, with information from
    // neighbouring face/cell.
    // Updates:
    //      - changedCell_, changedCells_
    //      - statistics: nEvals_, nUnvisitedCells_

    ++nEvals_;

    const bool wasValid = cellInfo.valid(td_);

    const bool propagate =
        cellInfo.updateCell
        (
            mesh_,
            celli,
            neighbourFacei,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        if (changedCell_.set(celli))
        {
            changedCells_.push_back(celli);
        }
    }

    if (!wasValid && cellInfo.valid(td_))
    {
        --nUnvisitedCells_;
    }

    return propagate;
}


template<class Type, class TrackingData>
bool Foam::FaceCellWave<Type, TrackingData>::updateFace
(
    const label facei,
    const label neighbourCelli,
    const Type& neighbourInfo,
    const scalar tol,
    Type& faceInfo
)
{
    // Update info for facei, at position pt, with information from
    // neighbouring face/cell.
    // Updates:
    //      - changedFace_, changedFaces_,
    //      - statistics: nEvals_, nUnvisitedFaces_

    ++nEvals_;

    const bool wasValid = faceInfo.valid(td_);

    const bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            facei,
            neighbourCelli,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        if (changedFace_.set(facei))
        {
            changedFaces_.push_back(facei);
        }
    }

    if (!wasValid && faceInfo.valid(td_))
    {
        --nUnvisitedFaces_;
    }

    return propagate;
}


template<class Type, class TrackingData>
bool Foam::FaceCellWave<Type, TrackingData>::updateFace
(
    const label facei,
    const Type& neighbourInfo,
    const scalar tol,
    Type& faceInfo
)
{
    // Update info for facei, at position pt, with information from
    // same face.
    // Updates:
    //      - changedFace_, changedFaces_,
    //      - statistics: nEvals_, nUnvisitedFaces_

    ++nEvals_;

    const bool wasValid = faceInfo.valid(td_);

    const bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            facei,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        if (changedFace_.set(facei))
        {
            changedFaces_.push_back(facei);
        }
    }

    if (!wasValid && faceInfo.valid(td_))
    {
        --nUnvisitedFaces_;
    }

    return propagate;
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::checkCyclic
(
    const polyPatch& patch
) const
{
    // For debugging: check status on both sides of cyclic

    const cyclicPolyPatch& nbrPatch =
        refCast<const cyclicPolyPatch>(patch).neighbPatch();

    forAll(patch, patchFacei)
    {
        const label i1 = patch.start() + patchFacei;
        const label i2 = nbrPatch.start() + patchFacei;

        if
        (
           !allFaceInfo_[i1].sameGeometry
            (
                mesh_,
                allFaceInfo_[i2],
                geomTol_,
                td_
            )
        )
        {
            FatalErrorInFunction
                << "   faceInfo:" << allFaceInfo_[i1]
                << "   otherfaceInfo:" << allFaceInfo_[i2]
                << abort(FatalError);
        }

        if (changedFace_.test(i1) != changedFace_.test(i2))
        {
            FatalErrorInFunction
                << "   faceInfo:" << allFaceInfo_[i1]
                << "   otherfaceInfo:" << allFaceInfo_[i2]
                << "   changedFace:" << changedFace_.test(i1)
                << "   otherchangedFace:" << changedFace_.test(i2)
                << abort(FatalError);
        }
    }
}


template<class Type, class TrackingData>
template<class PatchType>
bool Foam::FaceCellWave<Type, TrackingData>::hasPatch() const
{
    for (const polyPatch& p : mesh_.boundaryMesh())
    {
        if (isA<PatchType>(p))
        {
            return true;
        }
    }
    return false;
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::setFaceInfo
(
    const label facei,
    const Type& faceInfo
)
{
    const bool wasValid = allFaceInfo_[facei].valid(td_);

    // Copy info for facei
    allFaceInfo_[facei] = faceInfo;

    // Maintain count of unset faces
    if (!wasValid && allFaceInfo_[facei].valid(td_))
    {
        --nUnvisitedFaces_;
    }

    // Mark facei as visited and changed (both on list and on face itself)
    changedFace_.set(facei);
    changedFaces_.push_back(facei);
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::setFaceInfo
(
    const labelUList& changedFaces,
    const List<Type>& changedFacesInfo
)
{
    forAll(changedFaces, changedFacei)
    {
        const label facei = changedFaces[changedFacei];

        const bool wasValid = allFaceInfo_[facei].valid(td_);

        // Copy info for facei
        allFaceInfo_[facei] = changedFacesInfo[changedFacei];

        // Maintain count of unset faces
        if (!wasValid && allFaceInfo_[facei].valid(td_))
        {
            --nUnvisitedFaces_;
        }

        // Mark facei as changed, both on list and on face itself.
        changedFace_.set(facei);
        changedFaces_.push_back(facei);
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::mergeFaceInfo
(
    const polyPatch& patch,
    const label nFaces,
    const labelUList& changedFaces,
    const List<Type>& changedFacesInfo
)
{
    // Merge face information into member data

    for (label changedFacei = 0; changedFacei < nFaces; ++changedFacei)
    {
        const Type& newInfo = changedFacesInfo[changedFacei];
        const label patchFacei = changedFaces[changedFacei];

        const label meshFacei = patch.start() + patchFacei;

        Type& currInfo = allFaceInfo_[meshFacei];

        if (!currInfo.equal(newInfo, td_))
        {
            updateFace
            (
                meshFacei,
                newInfo,
                propagationTol_,
                currInfo
            );
        }
    }
}


template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::getChangedPatchFaces
(
    const polyPatch& patch,
    const label startFacei,
    const label nFaces,
    labelList& changedPatchFaces,
    List<Type>& changedPatchFacesInfo
) const
{
    // Construct compact patchFace change arrays for a (slice of a) single
    // patch. changedPatchFaces in local patch numbering.
    // Return length of arrays.
    label nChanged = 0;

    for (label i = 0; i < nFaces; ++i)
    {
        const label patchFacei = i + startFacei;
        const label meshFacei = patch.start() + patchFacei;

        if (changedFace_.test(meshFacei))
        {
            changedPatchFaces[nChanged] = patchFacei;
            changedPatchFacesInfo[nChanged] = allFaceInfo_[meshFacei];
            ++nChanged;
        }
    }

    return nChanged;
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::leaveDomain
(
    const polyPatch& patch,
    const label nFaces,
    const labelUList& faceLabels,
    List<Type>& faceInfo
) const
{
    // Handle leaving domain. Implementation referred to Type

    const vectorField& fc = mesh_.faceCentres();

    for (label i = 0; i < nFaces; ++i)
    {
        const label patchFacei = faceLabels[i];
        const label meshFacei = patch.start() + patchFacei;

        faceInfo[i].leaveDomain(mesh_, patch, patchFacei, fc[meshFacei], td_);
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::enterDomain
(
    const polyPatch& patch,
    const label nFaces,
    const labelUList& faceLabels,
    List<Type>& faceInfo
) const
{
    // Handle entering domain. Implementation referred to Type

    const vectorField& fc = mesh_.faceCentres();

    for (label i = 0; i < nFaces; ++i)
    {
        const label patchFacei = faceLabels[i];
        const label meshFacei = patch.start() + patchFacei;

        faceInfo[i].enterDomain(mesh_, patch, patchFacei, fc[meshFacei], td_);
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::transform
(
    const tensorField& rotTensor,
    const label nFaces,
    List<Type>& faceInfo
)
{
    // Transform. Implementation referred to Type

    if (rotTensor.size() == 1)
    {
        const tensor& T = rotTensor[0];

        for (label facei = 0; facei < nFaces; ++facei)
        {
            faceInfo[facei].transform(mesh_, T, td_);
        }
    }
    else
    {
        for (label facei = 0; facei < nFaces; ++facei)
        {
            faceInfo[facei].transform(mesh_, rotTensor[facei], td_);
        }
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::offset
(
    const polyPatch&,
    const label cycOffset,
    const label nFaces,
    labelList& faces
)
{
    // Offset mesh face.
    // Used for transferring from one cyclic half to the other.

    for (label facei = 0; facei < nFaces; ++facei)
    {
        faces[facei] += cycOffset;
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleProcPatches()
{
    // Transfer all the information to/from neighbouring processors

    const globalMeshData& pData = mesh_.globalData();

    // Which patches are processor patches
    const labelList& procPatches = pData.processorPatches();

    // Which processors this processor is connected to
    const labelList& neighbourProcs = pData.topology().procNeighbours();

    // Reduce communication by only sending non-zero data,
    // but with multiply-connected processor/processor
    // (eg, processorCyclic) also need to send zero information
    // to keep things synchronised

    // Reset buffers, initialize for registerSend() bookkeeping
    pBufs_.clear();
    pBufs_.initRegisterSend();


    // Information to send
    DynamicList<Type> sendFacesInfo;
    DynamicList<label> sendFaces;

    for (const label patchi : procPatches)
    {
        const auto& procPatch =
            refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchi]);

        const label nbrProci = procPatch.neighbProcNo();

        // Resize buffers
        sendFaces.resize_nocopy(procPatch.size());
        sendFacesInfo.resize_nocopy(procPatch.size());

        // Determine which faces changed on current patch
        const label nSendFaces = getChangedPatchFaces
        (
            procPatch,
            0,
            procPatch.size(),
            sendFaces,
            sendFacesInfo
        );

        // Shrink
        sendFaces.resize(nSendFaces);
        sendFacesInfo.resize(nSendFaces);

        // Adapt wallInfo for leaving domain
        leaveDomain
        (
            procPatch,
            nSendFaces,
            sendFaces,
            sendFacesInfo
        );

        // Send to neighbour
        {
            UOPstream toNbr(nbrProci, pBufs_);
            toNbr << sendFaces << sendFacesInfo;

            // Record if send is required (data are non-zero)
            pBufs_.registerSend(nbrProci, !sendFaces.empty());

            if (debug & 2)
            {
                Pout<< " Processor patch " << patchi << ' ' << procPatch.name()
                    << "  send:" << sendFaces.size() << " to proc:" << nbrProci
                    << endl;
            }
        }
    }

    // Limit exchange to involved procs
    // - automatically discards unnecessary (unregistered) sends
    pBufs_.finishedNeighbourSends(neighbourProcs);


    for (const label patchi : procPatches)
    {
        const auto& procPatch =
            refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchi]);

        const label nbrProci = procPatch.neighbProcNo();

        if (!pBufs_.recvDataCount(nbrProci))
        {
            // Nothing to receive
            continue;
        }


        labelList receiveFaces;
        List<Type> receiveFacesInfo;
        {
            UIPstream is(nbrProci, pBufs_);
            is >> receiveFaces >> receiveFacesInfo;
        }

        const label nReceiveFaces = receiveFaces.size();

        if (debug & 2)
        {
            Pout<< " Processor patch " << patchi << ' ' << procPatch.name()
                << "  recv:" << receiveFaces.size() << " from proci:"
                << nbrProci << endl;
        }

        // Apply transform to received data for non-parallel planes
        if (!procPatch.parallel())
        {
            transform
            (
                procPatch.forwardT(),
                nReceiveFaces,
                receiveFacesInfo
            );
        }

        // Adapt wallInfo for entering domain
        enterDomain
        (
            procPatch,
            nReceiveFaces,
            receiveFaces,
            receiveFacesInfo
        );

        // Merge received info
        mergeFaceInfo
        (
            procPatch,
            nReceiveFaces,
            receiveFaces,
            receiveFacesInfo
        );
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleCyclicPatches()
{
    // Transfer information across cyclic halves.

    for (const polyPatch& patch : mesh_.boundaryMesh())
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(patch);

        if (cpp)
        {
            const auto& cycPatch = *cpp;
            const auto& nbrPatch = cycPatch.neighbPatch();

            // Allocate buffers
            labelList receiveFaces(patch.size());
            List<Type> receiveFacesInfo(patch.size());

            // Determine which faces changed
            const label nReceiveFaces = getChangedPatchFaces
            (
                nbrPatch,
                0,
                nbrPatch.size(),
                receiveFaces,
                receiveFacesInfo
            );

            // Adapt wallInfo for leaving domain
            leaveDomain
            (
                nbrPatch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo
            );

            if (!cycPatch.parallel())
            {
                // Received data from other half
                transform
                (
                    cycPatch.forwardT(),
                    nReceiveFaces,
                    receiveFacesInfo
                );
            }

            if (debug & 2)
            {
                Pout<< " Cyclic patch "
                    << cycPatch.index() << ' ' << cycPatch.name()
                    << "  Changed : " << nReceiveFaces
                    << endl;
            }

            // Half2: Adapt wallInfo for entering domain
            enterDomain
            (
                cycPatch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo
            );

            // Merge into global storage
            mergeFaceInfo
            (
                cycPatch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo
            );

            if (debug)
            {
                checkCyclic(cycPatch);
            }
        }
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleAMICyclicPatches()
{
    // Transfer information across cyclicAMI halves.

    for (const polyPatch& patch : mesh_.boundaryMesh())
    {
        const cyclicAMIPolyPatch* cpp = isA<cyclicAMIPolyPatch>(patch);

        // Note:
        // - can either do owner and neighbour separately or have owner
        //   do both
        // - separately means that neighbour will receive updated owner
        //   properties which might be beneficial or involve extra work?

        if (cpp)
        {
            const auto& cycPatch = *cpp;
            const auto& nbrPatch = cycPatch.neighbPatch();

            List<Type> receiveInfo;

            {
                // Get nbrPatch data (so not just changed faces). Do not use
                // a slice here since the leaveDomain might change the values
                List<Type> sendInfo(nbrPatch.patchSlice(allFaceInfo_));

                if (!nbrPatch.parallel() || nbrPatch.separated())
                {
                    // Adapt sendInfo for leaving domain
                    const vectorField::subField fc = nbrPatch.faceCentres();
                    forAll(sendInfo, i)
                    {
                        sendInfo[i].leaveDomain(mesh_, nbrPatch, i, fc[i], td_);
                    }
                }

                // Transfer sendInfo to cycPatch
                combine<Type, TrackingData> cmb(*this, cycPatch);

                if (cycPatch.applyLowWeightCorrection())
                {
                    List<Type> defVals
                    (
                        cycPatch.patchInternalList(allCellInfo_)
                    );

                    cycPatch.interpolate(sendInfo, cmb, receiveInfo, defVals);
                }
                else
                {
                    cycPatch.interpolate(sendInfo, cmb, receiveInfo);
                }
            }

            // Apply transform to received data for non-parallel planes
            if (!cycPatch.parallel())
            {
                transform
                (
                    cycPatch.forwardT(),
                    receiveInfo.size(),
                    receiveInfo
                );
            }

            if (!cycPatch.parallel() || cycPatch.separated())
            {
                // Adapt receiveInfo for entering domain
                const vectorField::subField fc = cycPatch.faceCentres();
                forAll(receiveInfo, i)
                {
                    receiveInfo[i].enterDomain(mesh_, cycPatch, i, fc[i], td_);
                }
            }

            // Merge into global storage

            const auto areaFraction(patch.areaFraction());

            forAll(receiveInfo, i)
            {
                if (areaFraction && areaFraction()[i] <= 0.5)
                {
                    // not coupled
                    continue;
                }

                const label meshFacei = cycPatch.start()+i;

                const Type& newInfo = receiveInfo[i];

                Type& currInfo = allFaceInfo_[meshFacei];

                if (newInfo.valid(td_) && !currInfo.equal(newInfo, td_))
                {
                    updateFace
                    (
                        meshFacei,
                        newInfo,
                        propagationTol_,
                        currInfo
                    );
                }
            }
        }
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleExplicitConnections()
{
    changedBaffles_.clear();

    // Collect all/any changed information touching a baffle
    for (const labelPair& baffle : explicitConnections_)
    {
        const label f0 = baffle.first();
        const label f1 = baffle.second();

        if (changedFace_.test(f0))
        {
            // f0 changed. Update information on f1.
            changedBaffles_.push_back(taggedInfoType(f1, allFaceInfo_[f0]));
        }

        if (changedFace_.test(f1))
        {
            // f1 changed. Update information on f0.
            changedBaffles_.push_back(taggedInfoType(f0, allFaceInfo_[f1]));
        }
    }


    // Update other side with changed information

    for (const taggedInfoType& updated : changedBaffles_)
    {
        const label tgtFace = updated.first;
        const Type& newInfo = updated.second;

        Type& currInfo = allFaceInfo_[tgtFace];

        if (!currInfo.equal(newInfo, td_))
        {
            updateFace
            (
                tgtFace,
                newInfo,
                propagationTol_,
                currInfo
            );
        }
    }

    changedBaffles_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class TrackingData>
Foam::FaceCellWave<Type, TrackingData>::FaceCellWave
(
    const polyMesh& mesh,
    UList<Type>& allFaceInfo,
    UList<Type>& allCellInfo,
    TrackingData& td
)
:
    FaceCellWaveBase(mesh),

    explicitConnections_(),
    allFaceInfo_(allFaceInfo),
    allCellInfo_(allCellInfo),
    td_(td),
    changedBaffles_(2*explicitConnections_.size()),
    hasCyclicPatches_(hasPatch<cyclicPolyPatch>()),
    hasCyclicAMIPatches_
    (
        returnReduceOr(hasPatch<cyclicAMIPolyPatch>())
    ),
    nEvals_(0)
{
    if
    (
        allFaceInfo.size() != mesh_.nFaces()
     || allCellInfo.size() != mesh_.nCells()
    )
    {
        FatalErrorInFunction
            << "face and cell storage not the size of mesh faces, cells:" << nl
            << "    allFaceInfo   :" << allFaceInfo.size() << nl
            << "    mesh_.nFaces():" << mesh_.nFaces() << nl
            << "    allCellInfo   :" << allCellInfo.size() << nl
            << "    mesh_.nCells():" << mesh_.nCells() << endl
            << exit(FatalError);
    }
}


template<class Type, class TrackingData>
Foam::FaceCellWave<Type, TrackingData>::FaceCellWave
(
    const polyMesh& mesh,
    const labelUList& changedFaces,
    const List<Type>& changedFacesInfo,
    UList<Type>& allFaceInfo,
    UList<Type>& allCellInfo,
    const label maxIter,
    TrackingData& td
)
:
    FaceCellWaveBase(mesh),

    explicitConnections_(),
    allFaceInfo_(allFaceInfo),
    allCellInfo_(allCellInfo),
    td_(td),
    changedBaffles_(2*explicitConnections_.size()),
    hasCyclicPatches_(hasPatch<cyclicPolyPatch>()),
    hasCyclicAMIPatches_
    (
        returnReduceOr(hasPatch<cyclicAMIPolyPatch>())
    ),
    nEvals_(0)
{
    if
    (
        allFaceInfo.size() != mesh_.nFaces()
     || allCellInfo.size() != mesh_.nCells()
    )
    {
        FatalErrorInFunction
            << "face and cell storage not the size of mesh faces, cells:" << nl
            << "    allFaceInfo   :" << allFaceInfo.size() << nl
            << "    mesh_.nFaces():" << mesh_.nFaces() << nl
            << "    allCellInfo   :" << allCellInfo.size() << nl
            << "    mesh_.nCells():" << mesh_.nCells() << endl
            << exit(FatalError);
    }

    // Copy initial changed faces data
    setFaceInfo(changedFaces, changedFacesInfo);

    // Iterate until nothing changes
    const label iter = iterate(maxIter);

    if ((maxIter > 0) && (iter >= maxIter))
    {
        FatalErrorInFunction
            << "Maximum number of iterations reached. Increase maxIter." << nl
            << "    maxIter:" << maxIter << nl
            << "    nChangedCells:" << nChangedCells() << nl
            << "    nChangedFaces:" << nChangedFaces() << endl
            << exit(FatalError);
    }
}


template<class Type, class TrackingData>
Foam::FaceCellWave<Type, TrackingData>::FaceCellWave
(
    const polyMesh& mesh,
    const labelPairList& explicitConnections,
    const bool handleCyclicAMI,
    const labelUList& changedFaces,
    const List<Type>& changedFacesInfo,
    UList<Type>& allFaceInfo,
    UList<Type>& allCellInfo,
    const label maxIter,
    TrackingData& td
)
:
    FaceCellWaveBase(mesh),

    explicitConnections_(explicitConnections),
    allFaceInfo_(allFaceInfo),
    allCellInfo_(allCellInfo),
    td_(td),
    changedBaffles_(2*explicitConnections_.size()),
    hasCyclicPatches_(hasPatch<cyclicPolyPatch>()),
    hasCyclicAMIPatches_
    (
        handleCyclicAMI
     && returnReduceOr(hasPatch<cyclicAMIPolyPatch>())
    ),
    nEvals_(0)
{
    if
    (
        allFaceInfo.size() != mesh_.nFaces()
     || allCellInfo.size() != mesh_.nCells()
    )
    {
        FatalErrorInFunction
            << "face and cell storage not the size of mesh faces, cells:" << nl
            << "    allFaceInfo   :" << allFaceInfo.size() << nl
            << "    mesh_.nFaces():" << mesh_.nFaces() << nl
            << "    allCellInfo   :" << allCellInfo.size() << nl
            << "    mesh_.nCells():" << mesh_.nCells() << endl
            << exit(FatalError);
    }

    // Copy initial changed faces data
    setFaceInfo(changedFaces, changedFacesInfo);

    // Iterate until nothing changes
    const label iter = iterate(maxIter);

    if ((maxIter > 0) && (iter >= maxIter))
    {
        FatalErrorInFunction
            << "Maximum number of iterations reached. Increase maxIter." << nl
            << "    maxIter:" << maxIter << nl
            << "    nChangedCells:" << nChangedCells() << nl
            << "    nChangedFaces:" << nChangedFaces() << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::faceToCell()
{
    // Propagate face to cell

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();
    const label nInternalFaces = mesh_.nInternalFaces();

    for (const label facei : changedFaces_)
    {
        if (!changedFace_.test(facei))
        {
            FatalErrorInFunction
                << "Face " << facei
                << " not marked as having been changed"
                << abort(FatalError);
        }

        const Type& newInfo = allFaceInfo_[facei];

        // Evaluate all connected cells

        // Owner
        {
            const label celli = owner[facei];
            Type& currInfo = allCellInfo_[celli];

            if (!currInfo.equal(newInfo, td_))
            {
                updateCell
                (
                    celli,
                    facei,
                    newInfo,
                    propagationTol_,
                    currInfo
                );
            }
        }

        // Neighbour.
        if (facei < nInternalFaces)
        {
            const label celli = neighbour[facei];
            Type& currInfo = allCellInfo_[celli];

            if (!currInfo.equal(newInfo, td_))
            {
                updateCell
                (
                    celli,
                    facei,
                    newInfo,
                    propagationTol_,
                    currInfo
                );
            }
        }

        // Reset status of face
        changedFace_.unset(facei);
    }

    // Handled all changed faces by now
    changedFaces_.clear();

    if (debug & 2)
    {
        Pout<< " Changed cells            : " << nChangedCells() << endl;
    }

    // Number of changedCells over all procs
    return returnReduce(nChangedCells(), sumOp<label>());
}


template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::cellToFace()
{
    // Propagate cell to face

    const cellList& cells = mesh_.cells();

    for (const label celli : changedCells_)
    {
        if (!changedCell_.test(celli))
        {
            FatalErrorInFunction
                << "Cell " << celli << " not marked as having been changed"
                << abort(FatalError);
        }

        const Type& newInfo = allCellInfo_[celli];

        // Evaluate all connected faces

        const labelList& faceLabels = cells[celli];
        for (const label facei : faceLabels)
        {
            Type& currInfo = allFaceInfo_[facei];

            if (!currInfo.equal(newInfo, td_))
            {
                updateFace
                (
                    facei,
                    celli,
                    newInfo,
                    propagationTol_,
                    currInfo
                );
            }
        }

        // Reset status of cell
        changedCell_.unset(celli);
    }

    // Handled all changed cells by now
    changedCells_.clear();


    // Transfer across any explicitly provided internal connections
    handleExplicitConnections();

    if (hasCyclicPatches_)
    {
        handleCyclicPatches();
    }

    if (hasCyclicAMIPatches_)
    {
        handleAMICyclicPatches();
    }

    if (Pstream::parRun())
    {
        handleProcPatches();
    }

    if (debug & 2)
    {
        Pout<< " Changed faces            : " << nChangedFaces() << endl;
    }


    // Number of changedFaces over all procs
    return returnReduce(nChangedFaces(), sumOp<label>());
}


// Iterate
template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::iterate(const label maxIter)
{
    if (maxIter < 0)
    {
        return 0;
    }

    if (hasCyclicPatches_)
    {
        handleCyclicPatches();
    }

    if (hasCyclicAMIPatches_)
    {
        handleAMICyclicPatches();
    }

    if (Pstream::parRun())
    {
        handleProcPatches();
    }

    label iter = 0;

    for (/*nil*/; iter < maxIter; ++iter)
    {
        if (debug)
        {
            Info<< " Iteration " << iter << endl;
        }

        nEvals_ = 0;
        const label nCells = faceToCell();
        const label nFaces = nCells ? cellToFace() : 0;

        if (debug)
        {
            Info<< " Total evaluations     : "
                << nEvals_ << nl
                << " Changed cells / faces : "
                << nCells << " / " << nFaces << nl
                << " Pending cells / faces : "
                << nUnvisitedCells_ << " / " << nUnvisitedFaces_ << nl;
        }

        if (!nCells || !nFaces)
        {
            break;
        }
    }

    return iter;
}


// ************************************************************************* //
