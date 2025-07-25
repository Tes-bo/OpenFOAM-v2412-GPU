/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "globalIndexAndTransform.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::wallPointAddressing::setNull
(
    const globalIndexAndTransform& gt
)
{
    if (data_ == labelPair(-1, -1))
    {
        const label proci = 0;
        const label index = -1;

        // Add sending side transform
        data_ = gt.encode(proci, index, gt.nullTransformIndex());
    }
}


// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::wallPointAddressing::update
(
    const point& pt,
    const wallPointAddressing& w2,
    const scalar tol,
    TrackingData& td
)
{
    const scalar dist2 = magSqr(pt - w2.origin());

    if (valid(td))
    {
        const scalar diff = distSqr() - dist2;

        if (diff < 0)
        {
            // already nearer to pt
            return false;
        }

        if ((diff < SMALL) || ((distSqr() > SMALL) && (diff/distSqr() < tol)))
        {
            // don't propagate small changes
            return false;
        }
    }

    // Either *this is not yet valid or w2 is closer
    {
        // current not yet set so use any value
        distSqr() = dist2;
        origin() = w2.origin();
        data_ = w2.data();

        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wallPointAddressing::wallPointAddressing()
:
    wallPoint(),
    data_(labelPair(-1, -1))
{}


inline Foam::wallPointAddressing::wallPointAddressing
(
    const point& origin,
    const labelPair& data,
    const scalar distSqr
)
:
    wallPoint(origin, distSqr),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline void Foam::wallPointAddressing::leaveDomain
(
    const polyMesh& mesh,
    const polyPatch& pp,
    const label patchFacei,
    const point& faceCentre,
    TrackingData& td
)
{
    const auto& gt = mesh.globalData().globalTransforms();

    setNull(gt);

    // Convert to relative form
    wallPoint::leaveDomain(mesh, pp, patchFacei, faceCentre, td);

    // Add to transformation
    const label transformi = gt.transformIndex(data_);
    const label proci = gt.processor(data_);
    const label index = gt.index(data_);

    // Add sending side transform
    const label newTransformi = gt.addToTransformIndex
    (
        transformi,
        pp.index()
    );
    data_ = gt.encode(proci, index, newTransformi);
}


//// Update absolute geometric quantities. Note that distance (distSqr_)
//// is not affected by leaving/entering domain.
//template<class TrackingData>
//inline void Foam::wallPointAddressing::enterDomain
//(
//    const polyMesh& mesh,
//    const polyPatch& pp,
//    const label patchFacei,
//    const point& faceCentre,
//    TrackingData& td
//)
//{
//    // back to absolute form
//    wallPoint::enterDomain(mesh, pp, patchFacei, faceCentre, td);
//}


// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::wallPointAddressing::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label neighbourFacei,
    const wallPointAddressing& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    setNull(mesh.globalData().globalTransforms());

    const vectorField& cellCentres = mesh.primitiveMesh::cellCentres();

    return update
    (
        cellCentres[thisCelli],
        neighbourWallInfo,
        tol,
        td
    );
}


// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::wallPointAddressing::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label,
    const wallPointAddressing& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    setNull(mesh.globalData().globalTransforms());

    const vectorField& faceCentres = mesh.faceCentres();

    return update
    (
        faceCentres[thisFacei],
        neighbourWallInfo,
        tol,
        td
    );
}


// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::wallPointAddressing::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const wallPointAddressing& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    setNull(mesh.globalData().globalTransforms());

    const vectorField& faceCentres = mesh.faceCentres();

    return update
    (
        faceCentres[thisFacei],
        neighbourWallInfo,
        tol,
        td
    );
}


// ************************************************************************* //
