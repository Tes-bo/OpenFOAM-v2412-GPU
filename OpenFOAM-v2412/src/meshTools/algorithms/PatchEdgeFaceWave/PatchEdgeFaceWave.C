/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "PatchEdgeFaceWave.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update info for edgeI, at position pt, with information from
// neighbouring face.
// Updates:
//      - changedEdge_, changedEdges_,
//      - statistics: nEvals_, nUnvisitedEdges_
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData
>
bool Foam::PatchEdgeFaceWave<PrimitivePatchType, Type, TrackingData>::
updateEdge
(
    const label edgei,
    const label neighbourFacei,
    const Type& neighbourInfo,
    Type& edgeInfo
)
{
    nEvals_++;

    bool wasValid = edgeInfo.valid(td_);

    bool propagate =
        edgeInfo.updateEdge
        (
            mesh_,
            patch_,
            edgei,
            neighbourFacei,
            neighbourInfo,
            propagationTol_,
            td_
        );

    if (propagate)
    {
        if (changedEdge_.set(edgei))
        {
            changedEdges_.push_back(edgei);
        }
    }

    if (!wasValid && edgeInfo.valid(td_))
    {
        --nUnvisitedEdges_;
    }

    return propagate;
}


// Update info for facei, at position pt, with information from
// neighbouring edge.
// Updates:
//      - changedFace_, changedFaces_,
//      - statistics: nEvals_, nUnvisitedFace_
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData
>
bool Foam::PatchEdgeFaceWave<PrimitivePatchType, Type, TrackingData>::
updateFace
(
    const label facei,
    const label neighbourEdgeI,
    const Type& neighbourInfo,
    Type& faceInfo
)
{
    nEvals_++;

    bool wasValid = faceInfo.valid(td_);

    bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            patch_,
            facei,
            neighbourEdgeI,
            neighbourInfo,
            propagationTol_,
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


template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData
>
void Foam::PatchEdgeFaceWave<PrimitivePatchType, Type, TrackingData>::
syncEdges()
{
    const globalMeshData& globalData = mesh_.globalData();
    const mapDistribute& map = globalData.globalEdgeSlavesMap();
    const bitSet& cppOrientation = globalData.globalEdgeOrientation();

    // Convert patch-edge data into cpp-edge data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Construct with all data in consistent orientation
    List<Type> cppEdgeData(map.constructSize());

    forAll(patchEdges_, i)
    {
        label patchEdgeI = patchEdges_[i];
        label coupledEdgeI = coupledEdges_[i];

        if (changedEdge_.test(patchEdgeI))
        {
            const Type& data = allEdgeInfo_[patchEdgeI];

            // Patch-edge data needs to be converted into coupled-edge data
            // (optionally flipped) and consistent in orientation with
            // master of coupled edge (optionally flipped)
            bool sameOrientation =
            (
                sameEdgeOrientation_[i]
             == cppOrientation[coupledEdgeI]
            );

            cppEdgeData[coupledEdgeI].updateEdge
            (
                mesh_,
                patch_,
                data,
                sameOrientation,
                propagationTol_,
                td_
            );
        }
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
        updateOp<PrimitivePatchType, Type, TrackingData>
        (
            mesh_,
            patch_,
            propagationTol_,
            td_
        ),
        transformOp<PrimitivePatchType, Type, TrackingData>
        (
            mesh_,
            patch_,
            propagationTol_,
            td_
        )
    );


    // Back from cpp-edge to patch-edge data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(patchEdges_, i)
    {
        label patchEdgeI = patchEdges_[i];
        label coupledEdgeI = coupledEdges_[i];

        const Type& data = cppEdgeData[coupledEdgeI];

        if (data.valid(td_))
        {
            bool sameOrientation =
            (
                sameEdgeOrientation_[i]
             == cppOrientation[coupledEdgeI]
            );

            allEdgeInfo_[patchEdgeI].updateEdge
            (
                mesh_,
                patch_,
                data,
                sameOrientation,
                propagationTol_,
                td_
            );

            if (changedEdge_.set(patchEdgeI))
            {
                changedEdges_.push_back(patchEdgeI);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Iterate, propagating changedEdgesInfo across patch, until no change (or
// maxIter reached). Initial edge values specified.
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData
>
Foam::PatchEdgeFaceWave<PrimitivePatchType, Type, TrackingData>::
PatchEdgeFaceWave
(
    const polyMesh& mesh,
    const PrimitivePatchType& patch,
    const labelList& changedEdges,
    const List<Type>& changedEdgesInfo,

    UList<Type>& allEdgeInfo,
    UList<Type>& allFaceInfo,

    const label maxIter,
    TrackingData& td
)
:
    PatchEdgeFaceWaveBase(mesh, patch.nEdges(), patch.size()),
    patch_(patch),
    allEdgeInfo_(allEdgeInfo),
    allFaceInfo_(allFaceInfo),
    td_(td),
    nEvals_(0)
{
    // Calculate addressing between patch_ and mesh.globalData().coupledPatch()
    // for ease of synchronisation
    PatchTools::matchEdges
    (
        patch_,
        mesh_.globalData().coupledPatch(),

        patchEdges_,
        coupledEdges_,
        sameEdgeOrientation_
    );


    if (allEdgeInfo_.size() != patch_.nEdges())
    {
        FatalErrorInFunction
            << "size of edgeInfo work array is not equal to the number"
            << " of edges in the patch" << nl
            << "    edgeInfo   :" << allEdgeInfo_.size() << nl
            << "    patch.nEdges:" << patch_.nEdges() << endl
            << exit(FatalError);
    }
    if (allFaceInfo_.size() != patch_.size())
    {
        FatalErrorInFunction
            << "size of edgeInfo work array is not equal to the number"
            << " of faces in the patch" << nl
            << "    faceInfo   :" << allFaceInfo_.size() << nl
            << "    patch.size:" << patch_.size() << endl
            << exit(FatalError);
    }


    // Set from initial changed edges data
    setEdgeInfo(changedEdges, changedEdgesInfo);

    if (debug)
    {
        Pout<< "Seed edges                : " << nChangedEdges() << endl;
    }

    // Iterate until nothing changes
    label iter = iterate(maxIter);

    if ((maxIter > 0) && (iter >= maxIter))
    {
        FatalErrorInFunction
            << "Maximum number of iterations reached. Increase maxIter." << nl
            << "    maxIter:" << maxIter << nl
            << "    changedEdges:" << nChangedEdges() << nl
            << "    changedFaces:" << nChangedFaces() << endl
            << exit(FatalError);
    }
}


template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData
>
Foam::PatchEdgeFaceWave<PrimitivePatchType, Type, TrackingData>::
PatchEdgeFaceWave
(
    const polyMesh& mesh,
    const PrimitivePatchType& patch,
    UList<Type>& allEdgeInfo,
    UList<Type>& allFaceInfo,
    TrackingData& td
)
:
    PatchEdgeFaceWaveBase(mesh, patch.nEdges(), patch.size()),
    patch_(patch),
    allEdgeInfo_(allEdgeInfo),
    allFaceInfo_(allFaceInfo),
    td_(td),
    nEvals_(0)
{
    // Calculate addressing between patch_ and mesh.globalData().coupledPatch()
    // for ease of synchronisation
    PatchTools::matchEdges
    (
        patch_,
        mesh_.globalData().coupledPatch(),

        patchEdges_,
        coupledEdges_,
        sameEdgeOrientation_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Copy edge information into member data
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData
>
void Foam::PatchEdgeFaceWave<PrimitivePatchType, Type, TrackingData>::
setEdgeInfo
(
    const labelList& changedEdges,
    const List<Type>& changedEdgesInfo
)
{
    forAll(changedEdges, changedEdgeI)
    {
        label edgeI = changedEdges[changedEdgeI];

        bool wasValid = allEdgeInfo_[edgeI].valid(td_);

        // Copy info for edgeI
        allEdgeInfo_[edgeI] = changedEdgesInfo[changedEdgeI];

        // Maintain count of unset edges
        if (!wasValid && allEdgeInfo_[edgeI].valid(td_))
        {
            --nUnvisitedEdges_;
        }

        // Mark edgeI as changed, both on list and on edge itself.

        if (changedEdge_.set(edgeI))
        {
            changedEdges_.push_back(edgeI);
        }
    }
}


// Propagate information from face to edge. Return number of edges changed.
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData
>
Foam::label Foam::PatchEdgeFaceWave<PrimitivePatchType, Type, TrackingData>::
faceToEdge()
{
    changedEdges_.clear();
    changedEdge_ = false;

    for (const label facei : changedFaces_)
    {
        if (!changedFace_.test(facei))
        {
            FatalErrorInFunction
                << "face " << facei
                << " not marked as having been changed" << nl
                << "This might be caused by multiple occurrences of the same"
                << " seed edge." << abort(FatalError);
        }

        const Type& neighbourWallInfo = allFaceInfo_[facei];

        // Evaluate all connected edges
        const labelList& fEdges = patch_.faceEdges()[facei];

        forAll(fEdges, fEdgeI)
        {
            label edgeI = fEdges[fEdgeI];

            Type& currentWallInfo = allEdgeInfo_[edgeI];

            if (!currentWallInfo.equal(neighbourWallInfo, td_))
            {
                updateEdge
                (
                    edgeI,
                    facei,
                    neighbourWallInfo,
                    currentWallInfo
                );
            }
        }
    }


    syncEdges();


    if (debug)
    {
        Pout<< "Changed edges             : " << nChangedEdges() << endl;
    }

    return returnReduce(nChangedEdges(), sumOp<label>());
}


// Propagate information from edge to face. Return number of faces changed.
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData
>
Foam::label Foam::PatchEdgeFaceWave<PrimitivePatchType, Type, TrackingData>::
edgeToFace()
{
    changedFaces_.clear();
    changedFace_ = false;

    const labelListList& edgeFaces = patch_.edgeFaces();

    for (const label edgei : changedEdges_)
    {
        if (!changedEdge_.test(edgei))
        {
            FatalErrorInFunction
                << "edge " << edgei
                << " not marked as having been changed" << nl
                << "This might be caused by multiple occurrences of the same"
                << " seed edge." << abort(FatalError);
        }

        const Type& neighbourWallInfo = allEdgeInfo_[edgei];

        // Evaluate all connected faces

        for (const label facei : edgeFaces[edgei])
        {
            Type& currentWallInfo = allFaceInfo_[facei];

            if (!currentWallInfo.equal(neighbourWallInfo, td_))
            {
                updateFace
                (
                    facei,
                    edgei,
                    neighbourWallInfo,
                    currentWallInfo
                );
            }
        }
    }

    if (debug)
    {
        Pout<< "Changed faces             : " << nChangedFaces() << endl;
    }

    return returnReduce(nChangedFaces(), sumOp<label>());
}


// Iterate
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData
>
Foam::label Foam::PatchEdgeFaceWave<PrimitivePatchType, Type, TrackingData>::
iterate
(
    const label maxIter
)
{
    // Make sure coupled edges contain same info
    syncEdges();

    nEvals_ = 0;

    label iter = 0;

    while (iter < maxIter)
    {
        if (debug)
        {
            Pout<< "Iteration " << iter << endl;
        }

        label nFaces = edgeToFace();

        if (debug)
        {
            Pout<< "Total changed faces       : " << nFaces << endl;
        }

        if (nFaces == 0)
        {
            break;
        }

        label nEdges = faceToEdge();

        if (debug)
        {
            Pout<< "Total changed edges       : " << nEdges << nl
                << "Total evaluations         : " << nEvals_ << nl
                << "Remaining unvisited edges : " << nUnvisitedEdges_ << nl
                << "Remaining unvisited faces : " << nUnvisitedFaces_ << nl
                << endl;
        }

        if (nEdges == 0)
        {
            break;
        }

        iter++;
    }

    return iter;
}


// ************************************************************************* //
